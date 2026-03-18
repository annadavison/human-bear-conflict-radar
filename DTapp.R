library("leaflet")
library(leaflet.extras)
library(geosphere)
library(raster)
library(sf)
library(randomForest)

# general libraries
library(httr)
library(dplyr)
library(jsonlite)
library(shiny)
library(shinyjs)
library(lubridate)

# load the sensincluesr package
library(devtools)
devtools::install_github("sensingclues/sensingcluesr@v1.0.2", upgrade = "never")

# other necessities
BearIcon <- makeIcon("Bear1.png", iconWidth = 52, iconHeight = 65)
# Define the cutoff for the risk map
cutoff <- 0.15
# Define the ratio for new vs repeat offenders for the conflict risk map
ratio <- 0.5

# Load in the components for the dynamic conflict risk map
message("Loading the Bear HWC risk map")
conflictmodel <- readRDS("dynamic_conflict_model.rds")
predictDF <- read.csv("predict_riskmap.csv")
predmonthDF <- read.csv("month_prediction.csv")
mapbase <- raster("HWC risk map.tif")

# Load in the altitude map for the dispersal calculation
Altitude_map <- raster("DEM_Stara_Planina.tif")
Altitude_map <- Altitude_map*10 # dam to m

# Load the source library
source("bearDT-lib.R")
HotspotColor = c("grey80","#d9f0d3","#a6d96a","#ffff99","#fdae61","#e31a1c")

# Ensure that the risk is comparable
palette <- colorNumeric(
  palette = HotspotColor,
  domain = seq(0,1, length.out = 50),
  na.color = "transparent"
)


# login with functional user, credentials.json is required
credentials_path <- "/home/credentials_bear_radar_stara_planina.json"
if (file.exists(credentials_path)) {
  credentials <- jsonlite::fromJSON(credentials_path)
  cookie <- sensingcluesr::login_cluey(username = credentials$username,
                                       password = credentials$password)
} else {
  library(getPass)
  cookie <- sensingcluesr::login_cluey(getPass("Enter your Cluey username:"),
                                       getPass())}


# ----------------------------------------------------------------------------
# Parameters/variables
# ----------------------------------------------------------------------------

names <- c("Anoushka","Bernard","Christina","Djaner","Emiliya","Ferdinand",
           "Gergana","Hayes","Ivanka","John","Kathy","Lyudmil","Maya",
           "Nikola","Olga","Pavel","Queen","Rali","Stela","Tzvetan","Utha",
           "Vladimir","Wilma","Xanu","Yana","Zahari")
autoinvalidate <- 120               # number of seconds before recalculation starts
active_time <- 2                  # number of days before a bear becomes inactive

# environment
Env <- reactiveValues(
  time_now = Sys.time(),
  Bears = NULL,          # data frame with human bear conflict observations
  PAST.INCIDENTS = NULL, # data frame with all past incidents loaded into DT 
  risk_raster= NULL,     # raster layer with human bear conflict probabilities
  Bear.names = NULL,     # unique bear identifier
  Incident.IDs = NULL,   # unique incident identifier
  splines = NULL,        # spatial lines between bear and past incident(s)
  timepassed = NULL      # time since human bear conflict observation
)

# create initial map
map_start <- leaflet() %>% addProviderTiles("CartoDB.Positron") %>% flyToBounds(lng1 = 25.5, lng2 = 24, lat1 = 43.1, lat2 = 42.5) %>% 
  addLayersControl(overlayGroups = c("HWC risk map", "Active bears", "Past incidents","Current activity radius"), 
                   options = layersControlOptions(collapsed = FALSE), position = "topleft")

# ----------------------------------------------------------------------------
# Content
# ----------------------------------------------------------------------------
# time period to refresh data
autoInvalidate <- reactiveTimer(1000*autoinvalidate) # refresh every 120 seconds

observeEvent(autoInvalidate(),{
  
  # before loop
  time_now <- Sys.time()
  Bears <- Env$Bears
  PAST.INCIDENTS <- Env$PAST.INCIDENTS
  Bear.names <- Env$Bear.names
  splines <- NULL
  timepassed <- NULL
  
  ####################
  # get observations #
  ####################
  OBSERVATIONS <- NULL
  try(OBSERVATIONS <- getHWCobservations(cookie, from = (Sys.Date() - 2)), TRUE)
  if (is.null(OBSERVATIONS)) {message("FAILED GETTING OBSERVATIONS")} else {
    message("RUNNING CALCULATIONS")
    
    message("Removing old incidents")
    diftime <- difftime(time_now,OBSERVATIONS$datetime,units="days")
    OBSERVATIONS <- OBSERVATIONS[diftime<active_time,] # after 'active_time' days the bear location is too uncertain
    if (nrow(OBSERVATIONS)>0){
    if (is.null(PAST.INCIDENTS)) {
      firststep <- TRUE
      message("Running the first step")
      PAST.INCIDENTS <- data.frame(OBSERVATIONS)
      PAST.INCIDENTS$bear <- NA
    } else {
      firststep <- FALSE
      message("Checking for new incidents")
      #"check if OBSERVATIONS contain unique IDs not in PAST.INCIDENTS yet"
      check <- !OBSERVATIONS$entityId %in% PAST.INCIDENTS$entityId
      any_newincidents <- any(check)
      OBSERVATIONS <- OBSERVATIONS[check,]
      if(any_newincidents) {
        message("There are new incidents")
        NEW.INCIDENTS <- data.frame(OBSERVATIONS)
        NEW.INCIDENTS$bear <- NA
        PAST.INCIDENTS <- rbind(PAST.INCIDENTS,NEW.INCIDENTS)
      } else {
        message("No new incidents this time")
      }
    }
    }
    
    message("Fitting bears with incidents")
    
    if(is.null(Bears) & nrow(OBSERVATIONS)>0) {
      # Make the first bear:  
      Bears <- OBSERVATIONS[1,]
      if (firststep) {
        Bear.names <- names[1]
        Bears$name <- Bear.names
        PAST.INCIDENTS$bear[1] <- Bear.names
      } else {
        newnameindex <- max(which(names %in% Bear.names))+1
        troublemaker <- names[newnameindex]
        Bear.names <- c(Bear.names,troublemaker)
        Bears$name <- troublemaker
        Bear.latestInc <- Bears$entityId
        PAST.INCIDENTS$bear[PAST.INCIDENTS$entityId==Bear.latestInc] <- troublemaker
      }
      message("There is one bear in the dataset")
      message(paste("Its name is", last(Bear.names)))
      Bear.latestInc <- Bears$entityId
      
      # remove observation that belongs to the first bear:
      OBSERVATIONS <- OBSERVATIONS[OBSERVATIONS$entityId!=Bear.latestInc,]
    }
  }
  
  Ninc <- ifelse(is.null(OBSERVATIONS),0,nrow(OBSERVATIONS))
  if (Ninc > 0) {
    
    # There is more than one incident: loop through the incidents to check if there
    # are more bears:
    message("LOOPING THROUGH RECENT INCIDENTS") # make a function for this step
    for (i in 1:Ninc) {
      thisBear <- OBSERVATIONS[i,]
      troublemaker <- calcWhodunnit(incident=thisBear,
                                   suspects=Bears,
                                   Bear_names=Bears$name,
                                   Altitude_map=Altitude_map)
      new.bear <- troublemaker=="unknown"
      # add bear to df if there is a new bear related to a new incident
      if(new.bear) {
        newnameindex <- max(which(names %in% Bear.names))+1
        troublemaker <- names[newnameindex]
        Bear.names <- c(Bear.names,troublemaker)
        thisBear$name <- troublemaker
        Bears <- rbind(Bears,thisBear)
        message(paste("Bear",troublemaker, "is added to the bear dataset"))
      } else {
        thisBear$name <- troublemaker
        ### First remove latest sighting of troublemaker and then add new sighting of troublemaker
        Bears <- Bears[Bears$name!=troublemaker,]
        Bears <- rbind(Bears,thisBear)
        message(paste("There was a new incident with",troublemaker))
      }
      Bear.latestInc <- thisBear$entityId
      # add bear name to incident dataset
      PAST.INCIDENTS$bear[PAST.INCIDENTS$entityId==Bear.latestInc] <- troublemaker
      
    }
  } else {message("There are no new observations")}
  
  # remove inactive bears:
  if (!is.null(Bears)) {
    message("CHECKING FOR INACTIVE BEARS")
    timepassed <- as.numeric(difftime(time_now,Bears$datetime,unit="days"))
    do.remove <- timepassed>active_time
    if(any(do.remove)) {
      remove.names <- Bears$name[do.remove]
      Bears <- Bears[!do.remove,] # the bears that are not removed are left in the bear data set
      message(paste(paste(remove.names,"has been removed from the bear dataset due to inactivity."),collapse=" "))
    }
    
    if(nrow(Bears)>0) {
      message("GENERATING RISK MAP WITH LOCAL 48HR PREDICTIONS")
      # predict for 2 days ahead:
      BearRisk <- calcBearea(BearLoc = Bears,
                             time_now = time_now,
                             forecast = 48, # hours ahead
                             Altitude_map = Altitude_map
      )
      BearRisk <- projectRaster(BearRisk,mapbase,crs=crs(mapbase))
      month <- month(time_now)
      bearloc <- predmonthDF[,month]
      predictDF$bearloc <- bearloc
      predictDF <- cbind(predictDF,predict(conflictmodel, predictDF,"prob"))
      values(mapbase) <- predictDF$`TRUE`
      HWC_risk_map <- mapbase
      BearMap <- BearRisk * HWC_risk_map
      HWCrisk <- ((1-ratio)*BearMap) + (ratio*HWC_risk_map)
      HWCrisk[HWCrisk<cutoff] <- NA
      
      # make lines between active bears and past incidents
      lines <- Lines(lapply(Bears$name, function(x) {
        Line(PAST.INCIDENTS[PAST.INCIDENTS$bear==x, 3:4])}), ID="a")
      splines <- SpatialLines(list(lines))
      proj4string(splines) <- CRS("+init=epsg:4326")
      message("LINES MADE SUCCESSFULLY")
      
      # store passed time for the dispersal
      timepassed <- as.numeric(difftime(Sys.time(), Bears$datetime, unit="days"))
      
    } else {
      message("NO MORE ACTIVE BEARS; RISK MAP IS MONTHLY FORECAST")
      month <- month(time_now)
      bearloc <- predmonthDF[,month]
      predictDF$bearloc <- bearloc
      predictDF <- cbind(predictDF,predict(conflictmodel, predictDF,"prob"))
      values(mapbase) <- predictDF$`TRUE`
      HWCrisk <- ratio*mapbase
      HWCrisk[HWCrisk<cutoff] <- NA
    }
    
    if(nrow(Bears)==0) {Bears <- NULL}
  } else {
    message("NO RECENT INCIDENTS; RISK MAP IS MONTHLY FORECAST")
    month <- month(time_now)
    bearloc <- predmonthDF[,month]
    predictDF$bearloc <- bearloc
    predictDF <- cbind(predictDF,predict(conflictmodel, predictDF,"prob"))
    values(mapbase) <- predictDF$`TRUE`
    HWCrisk <- ratio*mapbase
    HWCrisk[HWCrisk<cutoff] <- NA
  }
  
  # after loop
  Env$time_now <- time_now
  Env$Bears <- Bears
  Env$PAST.INCIDENTS <- PAST.INCIDENTS
  Env$Bear.names <- Bear.names
  Env$risk_raster <- HWCrisk
  Env$splines <- splines
  Env$timepassed <- timepassed
  
})

# ----------------------------------------------------------------------------
# UI
# ----------------------------------------------------------------------------
# ui app
ui <- bootstrapPage(
  leafletOutput("MAP", height = '100vh'),
  absolutePanel(id = 'controls', top = 10, right = 10,
                wellPanel(
                  includeCSS("style.css"),
                  textOutput('connection'),
                  opacity = 0.75
                )
  ),
  useShinyjs(),
  # get timezone from browser
  tags$script("$(document).on('shiny:sessioninitialized', function(event) {
                                        var n = Intl.DateTimeFormat().resolvedOptions().timeZone;
                                        Shiny.onInputChange('user_timezone', n);});")
)

# ----------------------------------------------------------------------------
# Server
# ----------------------------------------------------------------------------
server <- function(input, output, session) {
  message("BEAR HWC APP SERVER START")
  
  # set up session data store
  session$userData$cookie <- NULL
  session$userData$authenticated <- FALSE
  session$userData$groups <- NULL
  
  # a modal for login to Focus
  dataModal <- function(failed = FALSE, failed_access = FALSE) {
    modalDialog(
      title = "Welcome to the Human Bear Conflict Radar",
      HTML("Please login with your Cluey credentials<br><br>"),
      textInput("username", "Cluey username"),
      passwordInput("password", "Cluey password"),
      size = "s",
      if (failed) 
        div(tags$b("Login failed, please try again", style = "color: red;")),
      if (failed_access) 
        div(tags$b("The bear radar has restricted access, if you would like more information please reach out to hello@sensingclues.org", style = "color: red;")),
      footer = tagList(
        actionButton("ok", "OK")
      )
    )
  }
  
  # This `observe` is suspended only with right user credential
  obs1 <- observe({
    showModal(dataModal())
  })
  
  # When OK button is pressed, attempt to login. If successful, remove the modal
  obs2 <- observe({
    req(input$ok)
    isolate({
      Username <- input$username
      Password <- input$password
    })
    
    session$userData$cookie <- sensingcluesr::login_cluey(username = Username, password = Password)
    if (!is.null(isolate(session$userData$cookie))) {
      session$userData$groups <- sensingcluesr::get_groups(session$userData$cookie)
      # the bear radar requires the cluey user to be part of the bear-radar-demo group
      if ("bear-radar-demo" %in% session$userData$groups$value) {
        obs1$suspend()
        removeModal()
        isolate(session$userData$authenticated <- TRUE)
      } else {
        # inform user that access is restricted
        showModal(dataModal(failed_access = TRUE))
      }
    } else {
      # inform user that login failed
      showModal(dataModal(failed = TRUE))
    }     
  })
  
  # information for user
  observe({
    message(paste("Timezone of the user is", input$user_timezone))
    if (!is.null(input$user_timezone)) {
      output$connection <- renderText(paste('Last prediction was at',
                                            format(Env$time_now, 
                                                   format = "%Y-%m-%d %H:%M:%S %Z",
                                                   tz = input$user_timezone)))
    }
  })
  
  observeEvent(Env$risk_raster, {
    
    if (!is.null(Env$risk_raster) & session$userData$authenticated) {
      message("RETRIEVEING THE MAP FROM SERVER")
      
      Bears <- Env$Bears
      HWCrisk <- Env$risk_raster
      PAST.INCIDENTS <- Env$PAST.INCIDENTS
      splines <- Env$splines
      timepassed <- Env$timepassed
      # hourly:
      timepassed <- timepassed*24
      if ((!is.null(Bears)) & # is there bear data at all?
          ifelse(is.null(Bears), FALSE, nrow(Bears) > 0) # have all the bears been removed because of inactivity?
      ) {
      Altitude <- extract(Altitude_map,Bears)
      
      upperrange_poly <- calcActivityArea(Bears[timepassed>=1,],dtime=timepassed[timepassed>=1],
                                          BearAlt=Altitude[timepassed>=1],
                                          bounds="upper",
                                          Altitude_map = Altitude_map)
      prefrange_poly <- calcActivityArea(Bears[timepassed>=1,],dtime=timepassed[timepassed>=1],
                                         BearAlt=Altitude[timepassed>=1],
                                         Altitude_map = Altitude_map)
      # ADDED THE CONVERSION TO SF TO BE ABLE TO PLOT
      upperrange_poly <- st_as_sf(upperrange_poly)
      prefrange_poly <- st_as_sf(prefrange_poly)
      b <- 1/0.04721162
      a <- 0.7444089*b
      alpha <- 2.0892
      beta <- 0.4872
      distanceU <- a*timepassed/(b+timepassed)
      distanceU[timepassed<24] <- exp(alpha+beta*log(timepassed/24))[timepassed<24] # use exponential dispersal model in first 24 hours
      distanceU <- distanceU*1000
      
      b <- 1/0.04161593
      a <- 0.2527037*b
      alpha <- 1.08046
      beta <- 0.48566
      distanceM <- a*timepassed/(b+timepassed)
      distanceM[timepassed<24] <- exp(alpha+beta*log(timepassed/24))[timepassed<24] # use exponential dispersal model in first 24 hours
      distanceM <- distanceM*1000
        
        leafletProxy("MAP", session) %>%
          addRasterImage(x=HWCrisk, 
                         col=palette,
                         opacity = .7,
                         layerId = "riskmap",group = "HWC risk map") %>% 
          addLegend("bottomleft", colors = c("#d9f0d3","#a6d96a","#ffff99","#fdae61","#e31a1c"),
                    labels = c("Negligible risk", "Very small risk", "Small risk","Medium risk","High risk"),
                    layerId = "legend", title = "Conflict risk levels", opacity = 0.8) %>% 
          addPolygons(data=upperrange_poly,color = "blue",
                      layerId = "Radius_LO",
                      group = "Current activity radius") %>% 
          addPolygons(data=prefrange_poly,color="red",
                      layerId = "Radius_HI",
                      group = "Current activity radius") %>% 
          addPolylines(data = splines, layerId = "lines", weight = 3,
                       group = "Past incidents")%>% 
          addMarkers(data=Bears, 
                     popup=paste0(
                       # "<strong>Bear name: </strong>", Bears$name, "<br>",
                       "<strong>Last incident: </strong>", Bears$datetime #, "<br>",
                       #"<strong>Tags: </strong>", HWCdetails$tags[Incident.IDs]
                     ),
                     layerId = paste("Bear",Bears$name),
                     group = "Active bears",
                     icon=BearIcon) %>% 
          addCircles(data=Bears,layerId = paste("Radius_LO",Bears$name),
                     group = "Current activity radius",
                     radius = distanceU,
                     fill = NA, 
                     weight = 2) %>% 
          addCircles(data=Bears,layerId = paste("Radius_HI",Bears$name),
                     group = "Current activity radius",
                     radius = distanceM,
                     fill = NA, color="red", 
                     weight = 2)
      } else {
        leafletProxy("MAP", session) %>%
          addRasterImage(x=HWCrisk, 
                         col=palette,
                         opacity = .7,
                         layerId = "riskmap",group = "HWC risk map") %>%
          addLegend("bottomleft", colors = c("#d9f0d3","#a6d96a","#ffff99","#fdae61","#e31a1c"),
                    labels = c("Negligible risk", "Very small risk", "Small risk","Medium risk","High risk"),
                    layerId = "legend", title = "Conflict risk levels", opacity = 0.8)
      }
    }
    message("FINISHED")
    
  })
  
  output$MAP <- renderLeaflet({
    map_start
  })
}

shinyApp(ui, server)

########