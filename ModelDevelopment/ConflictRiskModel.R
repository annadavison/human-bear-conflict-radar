##########################################
## K-FOLD CROSS VALIDATION FOR RISK MAP ##
##########################################

library(getPass)
library("sp")
library("leaflet")
library(leaflet.extras)
library("sensingcluesr")
library(geosphere)
library(raster)
library(httr)
library(dplyr)
library(shiny)
library(shinyjs)
library(lubridate)
library(loggit)

set.seed(20120)

HotspotColor = c("grey80","#d9f0d3","#a6d96a","#ffff99","#fdae61","#e31a1c")

palette <- colorNumeric(
  palette = HotspotColor,
  domain = seq(0,1, length.out = 50),
  na.color = "transparent"
)

#-------------------------------------------------------------------------------
#_______________________________________________________________________________
# Load in the Points from Sensing Clues
#-------------------------------------------------------------------------------
#_______________________________________________________________________________

cookie <- login_cluey(getPass("username"),
                      getPass("password"))
groups <- get_groups(cookie)
head(groups)

#-------------------------------------------------------------------------------

# Import the Signs of Bear Presence
# The concepts are bear or animal observation
# These points have longitude, latitude, time and species.
df2 <- get_observations(cookie, 
                        group = c("focus-project-8480276"),
                        from = "2000-01-01", to = "2024-12-31", allAttributes = TRUE,
                        filteredConcepts = c("https://sensingclues.poolparty.biz/SCCSSOntology/186","https://sensingclues.poolparty.biz/SCCSSOntology/266"))
SOPpoints <- df2 %>% dplyr::filter(grepl("bear", tags)) %>% 
  dplyr::filter(!duplicated(entityId)) %>% 
  dplyr::mutate(datetime = lubridate::as_datetime(when, tz = "UTC"), # format="%Y-%m-%dT%H:%M:%OS"
                lon = longitude,
                lat = latitude) %>% 
  dplyr::arrange(datetime) %>% 
  dplyr::select(any_of(c("lon","lat","entityId","datetime","longitude","latitude","tags", "where")))

N <- nrow(SOPpoints)
for(i in 1:N) {
  coords <- jsonlite::fromJSON(SOPpoints$where[i])$coordinates
  SOPpoints$lon[i] <- coords[1,1]
  SOPpoints$lat[i] <- coords[2,1]
}
coordinates(SOPpoints) <- ~lon+lat
proj4string(SOPpoints) <- CRS("+init=epsg:4326")
# Make a map of the bear observations
leaflet() %>% addTiles() %>% 
  #addHeatmap(data=track_coordinates,max = 25, radius = 20, blur = 20) %>% 
  addMarkers(data=SOPpoints, popup=SOPpoints$tags) %>% 
  addMeasure(primaryLengthUnit = "meters",primaryAreaUnit = "sqmeters")

SOPpoints$datetime
SOPpoints <- SOPpoints[year(SOPpoints$datetime) < 2025,]

#-------------------------------------------------------------------------------
# human-bear-conflict
# Download all the observations (points) which come under HWC or bears
# This data has information on longitude, latitude, time, species and description of conflict type
HWC <- get_observations(cookie, 
                        group = c("focus-project-8480276"),
                        from = "2000-01-01", to = "2024-12-31", allAttributes = TRUE,
                        filteredConcepts = c("https://sensingclues.poolparty.biz/SCCSSOntology/106","https://sensingclues.poolparty.biz/SCCSSOntology/266"))

HWC <- HWC %>% dplyr::filter(grepl("bear", tags)) %>% 
  dplyr::filter(!duplicated(entityId)) %>% 
  dplyr::mutate(datetime = lubridate::as_datetime(when, tz = "UTC"), # format="%Y-%m-%dT%H:%M:%OS"
                lon = longitude,
                lat = latitude) %>% 
  dplyr::arrange(datetime) %>% 
  dplyr::select(any_of(c("lon","lat","entityId","datetime","longitude","latitude","tags", "where", "projectName")))

N <- nrow(HWC)
# Add lon and lat to the HWC data frame
for(i in 1:N) {
  coords <- jsonlite::fromJSON(HWC$where[i])$coordinates
  HWC$lon[i] <- coords[1,1]
  HWC$lat[i] <- coords[2,1]
}
# Make the dataframe into coordinates I believe using the sp package
coordinates(HWC) <- ~lon+lat
# Set the coordinate system
proj4string(HWC) <- CRS("+init=epsg:4326")

# Plot the HWC points
leaflet() %>% addTiles() %>% 
  #addHeatmap(data=track_coordinates,max = 25, radius = 20, blur = 20) %>% 
  addMarkers(data=HWC, popup=HWC$tags) %>% 
  addMeasure(primaryLengthUnit = "meters",primaryAreaUnit = "sqmeters")
# Take a look at when the points were made
HWC$datetime
# Remove after 2024 (I kind of already did that but doesn't hurt to remove here too)
HWC <- HWC[year(HWC$datetime) < 2025,]


### Now take a look at the official data points ###
HWC_upload <- read.csv("CB_Bear_damages.csv", sep = ";")

# Make the dataframe into coordinates I believe using the sp package
coordinates(HWC_upload) <- ~Long+Lat
# Set the coordinate system
proj4string(HWC_upload) <- CRS("+init=epsg:4326")
# Plot the HWC points
leaflet() %>% addTiles() %>% 
  #addHeatmap(data=track_coordinates,max = 25, radius = 20, blur = 20) %>% 
  addMarkers(data=HWC_upload) %>% 
  addMeasure(primaryLengthUnit = "meters",primaryAreaUnit = "sqmeters")
# Create a as.POSIXlt column with the time information
HWC_upload@data$datetime <- as.POSIXct(HWC_upload@data$Date,format="%d.%m.%Y",tz="EET")

# Make combined HWC points to use in the risk map generation
# Keep just the datetime and tags for each for merging
HWC@data$Source <- rep("App", nrow(HWC))
HWC_n <- HWC[,c("datetime", "tags", "Source")]
names(HWC_upload@data)[names(HWC_upload@data) == "Type"] <- "tags"
HWC_up <- HWC_upload[,c("datetime", "tags", "Source")]
HWCpoints <- rbind(HWC_up, HWC_n)

#--------------------
## Load in the HSM ##
#--------------------

# The HSM was created using the methodology described in de Koning et al., (in review)
# which will be published in MethodsX and linked to the HBC Radar paper where further
# methods are described for this specific study
filename <- "HSM_gentoFeb2025.tif"
HSM <- raster(filename)

# Make circ filter function

# This is a custom function to smooth out the edges of the prior map in order
# to display it with less hard boundaries.
make_circ_filter<-function(radius, res){
  circ_filter<-matrix(NA, nrow=1+(2*radius/res), ncol=1+(2*radius/res))
  dimnames(circ_filter)[[1]]<-seq(-radius, radius, by=res)
  dimnames(circ_filter)[[2]]<-seq(-radius, radius, by=res)
  sweeper<-function(mat){
    for(row in 1:nrow(mat)){
      for(col in 1:ncol(mat)){
        dist<-sqrt((as.numeric(dimnames(mat)[[1]])[row])^2 +
                     (as.numeric(dimnames(mat)[[1]])[col])^2)
        if(dist<=radius) {mat[row, col]<-1}
      }
    }
    return(mat)
  }
  out<-sweeper(circ_filter)
  return(out)
}

# Cut the HWC points to the HSM
transf <- HWCpoints
# Reproject raster to Decimal Degrees
RP_HSM <- projectRaster(HSM, crs = CRS("+proj=longlat +datum=WGS84"))
extent <- extent(RP_HSM)
crs(transf) <- crs(RP_HSM)
clip <- crop(transf, RP_HSM, inverse = F)
leaflet() %>% addTiles() %>% addMarkers(data = clip)

# These are the final HWC points to use:
HWCpoints <- clip


####################################
## determine landscape patchiness ##
####################################

library(landscapemetrics)
library(lidaRtRee)
# make raster of very suitable, somewhat suitable and unsuitable habitat:
classification <- (HSM > 10) * 2 + (HSM < 10 & HSM > 5)
edges <- get_boundaries(classification)[[1]]
plot(edges)
edges[is.na(edges)] <- 0
filterradius2 <- 2000 # is the daily displacement of 50% of bears
cf2<-make_circ_filter(filterradius2, 100)
edges_persqkm <-focal(edges[[1]], w=cf2, fun=function(x){mean(x, na.rm=TRUE)})
plot(edges_persqkm)
edges_persqkm <- convert_raster(edges_persqkm, pkg = "raster")

####################################################
## determine core area and distance to core area: ##
####################################################

# Now the core area determined within the workshop is now used instead of the
# one generated purely from the HSM. The workshop version removed irrelevant 
# sections and expanded one of them
filename <- "bearcore.tif"
core_area <- raster::raster(filename)
dist_to_core <- distance(core_area)

########################
## Temporal variables ##
########################
# First method utilises percentages of activity within and outside the core bear
# area derived from bear GPS data. This data was aggregated to a daily mean location
# for each bear, this was then grouped by month for all bears and percentages
# within and outside the core areas determined and listed in a CSV file.

# Load in the percentages of bears in and out of the core_area
monthper <- read.csv("month_percentage.csv")

# Then load in the other method for the temporal variables which uses the 
# activity within conflict prone zones

# This was determined in a similar way to the previous method instead this time the
# daily aggregated mean location was determined to be either within predetermined 
# conflict prone habitats (CLC 242 & 243) or outside. A percentage of total points
# within the conflict prone zones is taken for each month and the areas outside these
# zones are assigned a value of 0.
monthact <- read.csv("MonthlyActivityConflictZones_reshuff.csv", sep = ";")
# Load in the raster mask for the conflict zones
filename <- "dangerzone_242_243.tif"
confzone <- raster(filename)

################################################################################
#-------------------------------------------------------------------------------
# COMPARISON SETUP #
#-------------------------------------------------------------------------------
################################################################################

# So there are a few elements to compare for this selection of the perfect 
# conflict risk map

# Here is the layout of the experiments for clarity:

# IN/OUT CORE ACTIVITY as temporal input
###### Pseudeo-absence 3km BUFFER
########### Unconstrained
########### 25 Node cap
###### Pseudo-absence RANDOM
########### Unconstrained
########### 25 Node cap
###### Pseudo-absence 5km BUFFER
########### Unconstrained
########### 25 Node cap

# Train with a Random forest, logistic regression and GBM

# Then test with K-fold cross validation and spatial blocking CV

# This is repeated, first with the alternative variable CONFLICT ZONE ACTIVITY
# Then with the DISTANCE FROM CORE BEAR AREA REMOVED


################################################################################
# IN/OUT CORE ACTIVITY as temporal input
################################################################################
################################################################################
# IN/OUT CORE ACTIVITY as temporal input
################################################################################
# Reset from last
HWCpoints <- clip
######################
## Assign Variables ##
######################

# Habitat Suitability
HWCpoints$habitat_suitability <- raster::extract(HSM,HWCpoints)

# Edges per square km
HWCpoints$edges_persqkm <- raster::extract(edges_persqkm,HWCpoints)
hist(edges_persqkm[edges_persqkm>0],xlim=c(0,0.5))
hist(HWCpoints$edges_persqkm,xlim=c(0,0.5), breaks=10)

# Distance to the Core
HWCpoints$dist_to_core <- raster::extract(dist_to_core,HWCpoints)
hist(dist_to_core, xlim=c(0,60000))
hist(HWCpoints$dist_to_core, xlim=c(0,60000),breaks=15)

# Monthly activity within and outside the core
HWCpoints@data$bearloc <- raster::extract(core_area, HWCpoints)
for (month in c(1:12)){
  HWCpoints@data$bearloc[is.na(HWCpoints@data$bearloc) & (month(HWCpoints@data$datetime) == month)] <- monthper$Percentage[month]
  HWCpoints@data$bearloc[(HWCpoints@data$bearloc == 1) & (month(HWCpoints@data$datetime) == month)] <- 100- monthper$Percentage[month]
}

# Monthly activity for the conflict prone zones
monthact <- read.csv("MonthlyActivityConflictZones_reshuff.csv", sep = ";")
# Load in the raster mask for the conflict zones
filename <- "dangerzone_242_243.tif"
confzone <- raster(filename)

# Since there are 0 values for this new method and we need 0 to also be outside
# we introduce a binary raster of inside or outside the conflict zone.
HWCpoints@data$bearloc2 <- raster::extract(confzone, HWCpoints)
HWCpoints@data$bearact <- rep(NA, nrow(HWCpoints))
for (month in c(1:12)){
  HWCpoints@data$bearact[is.na(HWCpoints@data$bearloc2) & (month(HWCpoints@data$datetime) == month)] <- 0
  HWCpoints@data$bearact[(HWCpoints@data$bearloc2 == 1) & (month(HWCpoints@data$datetime) == month)] <- monthact$Activity[month]
}
HWCpoints@data$bearloc2[is.na(HWCpoints@data$bearloc2)] <- 0


######################################
## Generate the pseudo-absence data ##
######################################

#-------------------------
# Met1 - Original Method #
#-------------------------
buffer <- buffer(HWCpoints, width=3000) # draw 3km buffer around HWC points
plot(buffer)
plot(SOPpoints,add=TRUE)
plot(SOPpoints[buffer,],add=TRUE,col="blue")
plot(SOPpoints[is.na(over(SOPpoints,buffer)),],add=TRUE,col="red")
# generate pseudo absence data:
# Filters down to roughly the same number of points as the HWCpoints
set.seed(20120)
met1_PABpoints <- data.frame(coordinates(core_area)[sample(1:ncell(core_area),200),])
coordinates(met1_PABpoints) <- ~x+y
proj4string(met1_PABpoints) <- proj4string(core_area)
# extract features:
met1_PABpoints$dist_to_core <- raster::extract(dist_to_core,met1_PABpoints)
met1_PABpoints$edges_persqkm <- raster::extract(edges_persqkm,met1_PABpoints)
met1_PABpoints$habitat_suitability <- raster::extract(HSM,met1_PABpoints)
met1_PABpoints <- spTransform(met1_PABpoints,CRS("+init=epsg:4326"))
met1_PABpoints <- met1_PABpoints[is.na(over(met1_PABpoints,buffer)),] # select only points that fall outside buffer area
# Add the month and then the corresponding percentages
# CURRENTLY JUST RANDOM
set.seed(20120)
met1_PABpoints$month <- sample(c(1:12), nrow(met1_PABpoints), replace = TRUE)
met1_PABpoints@data$bearloc <- raster::extract(core_area, met1_PABpoints)
for (month in c(1:12)){
  met1_PABpoints@data$bearloc[is.na(met1_PABpoints@data$bearloc) & (met1_PABpoints@data$month == month)] <- monthper$Percentage[month]
  met1_PABpoints@data$bearloc[(met1_PABpoints@data$bearloc == 1) & (met1_PABpoints@data$month == month)] <- 100-monthper$Percentage[month]
}
# Now for the alternative variable, confzone
met1_PABpoints@data$bearloc2 <- raster::extract(confzone, met1_PABpoints)
met1_PABpoints@data$bearact <- rep(NA, nrow(met1_PABpoints))
for (month in c(1:12)){
  met1_PABpoints@data$bearact[is.na(met1_PABpoints@data$bearloc2) & (met1_PABpoints@data$month == month)] <- 0
  met1_PABpoints@data$bearact[(met1_PABpoints@data$bearloc2 == 1) & (met1_PABpoints@data$month == month)] <- monthact$Activity[month]
}
met1_PABpoints@data$bearloc2[is.na(met1_PABpoints@data$bearloc2)] <- 0

leaflet() %>% addTiles() %>%
  addMarkers(data=met1_PABpoints) %>% 
  addMeasure(primaryLengthUnit = "meters",primaryAreaUnit = "sqmeters")


#---------------------------------
# Met2 - Random generation #
#---------------------------------

# No distance-based filtering process, randomly generated

# generate pseudo absence data:
# Filters down to roughly the same number of points as the HWCpoints
set.seed(20120)
met2_PABpoints <- data.frame(coordinates(core_area)[sample(1:ncell(core_area),200),])
coordinates(met2_PABpoints) <- ~x+y
proj4string(met2_PABpoints) <- proj4string(core_area)
# extract features:
met2_PABpoints$dist_to_core <- raster::extract(dist_to_core,met2_PABpoints)
met2_PABpoints$edges_persqkm <- raster::extract(edges_persqkm,met2_PABpoints)
met2_PABpoints$habitat_suitability <- raster::extract(HSM,met2_PABpoints)
met2_PABpoints <- spTransform(met2_PABpoints,CRS("+init=epsg:4326"))
# Add the month and then the corresponding percentages
# CURRENTLY JUST RANDOM
set.seed(20120)
met2_PABpoints$month <- sample(c(1:12), nrow(met2_PABpoints), replace = TRUE)
met2_PABpoints@data$bearloc <- raster::extract(core_area, met2_PABpoints)
for (month in c(1:12)){
  met2_PABpoints@data$bearloc[is.na(met2_PABpoints@data$bearloc) & (met2_PABpoints@data$month == month)] <- monthper$Percentage[month]
  met2_PABpoints@data$bearloc[(met2_PABpoints@data$bearloc == 1) & (met2_PABpoints@data$month == month)] <- 100-monthper$Percentage[month]
}
# Now for the alternative variable, confzone
met2_PABpoints@data$bearloc2 <- raster::extract(confzone, met2_PABpoints)
met2_PABpoints@data$bearact <- rep(NA, nrow(met2_PABpoints))
for (month in c(1:12)){
  met2_PABpoints@data$bearact[is.na(met2_PABpoints@data$bearloc2) & (met2_PABpoints@data$month == month)] <- 0
  met2_PABpoints@data$bearact[(met2_PABpoints@data$bearloc2 == 1) & (met2_PABpoints@data$month == month)] <- monthact$Activity[month]
}
met2_PABpoints@data$bearloc2[is.na(met2_PABpoints@data$bearloc2)] <- 0



leaflet() %>% addTiles() %>%
  addMarkers(data=met2_PABpoints) %>% 
  addMeasure(primaryLengthUnit = "meters",primaryAreaUnit = "sqmeters")

#----------------------------
# Met3 - Increase the buffer
#----------------------------
# Only difference is that the buffer is increased to 5km from 3km

buffer <- buffer(HWCpoints, width=5000) # draw 5km buffer around HWC points
plot(buffer)
plot(SOPpoints,add=TRUE)
plot(SOPpoints[buffer,],add=TRUE,col="blue")
plot(SOPpoints[is.na(over(SOPpoints,buffer)),],add=TRUE,col="red")
# generate pseudo absence data:
# Filters down to roughly the same number of points as the HWCpoints
set.seed(20120)
met3_PABpoints <- data.frame(coordinates(core_area)[sample(1:ncell(core_area),200),])
coordinates(met3_PABpoints) <- ~x+y
proj4string(met3_PABpoints) <- proj4string(core_area)
# extract features:
met3_PABpoints$dist_to_core <- raster::extract(dist_to_core,met3_PABpoints)
met3_PABpoints$edges_persqkm <- raster::extract(edges_persqkm,met3_PABpoints)
met3_PABpoints$habitat_suitability <- raster::extract(HSM,met3_PABpoints)
met3_PABpoints <- spTransform(met3_PABpoints,CRS("+init=epsg:4326"))
met3_PABpoints <- met3_PABpoints[is.na(over(met3_PABpoints,buffer)),] # select only points that fall outside buffer area
# Add the month and then the corresponding percentages
# CURRENTLY JUST RANDOM
set.seed(20120)
met3_PABpoints$month <- sample(c(1:12), nrow(met3_PABpoints), replace = TRUE)
met3_PABpoints@data$bearloc <- raster::extract(core_area, met3_PABpoints)
for (month in c(1:12)){
  met3_PABpoints@data$bearloc[is.na(met3_PABpoints@data$bearloc) & (met3_PABpoints@data$month == month)] <- monthper$Percentage[month]
  met3_PABpoints@data$bearloc[(met3_PABpoints@data$bearloc == 1) & (met3_PABpoints@data$month == month)] <- 100-monthper$Percentage[month]
}
# Now for the alternative variable, confzone
met3_PABpoints@data$bearloc2 <- raster::extract(confzone, met3_PABpoints)
met3_PABpoints@data$bearact <- rep(NA, nrow(met3_PABpoints))
for (month in c(1:12)){
  met3_PABpoints@data$bearact[is.na(met3_PABpoints@data$bearloc2) & (met3_PABpoints@data$month == month)] <- 0
  met3_PABpoints@data$bearact[(met3_PABpoints@data$bearloc2 == 1) & (met3_PABpoints@data$month == month)] <- monthact$Activity[month]
}
met3_PABpoints@data$bearloc2[is.na(met3_PABpoints@data$bearloc2)] <- 0

leaflet() %>% addTiles() %>%
  addMarkers(data=met3_PABpoints) %>% 
  addMeasure(primaryLengthUnit = "meters",primaryAreaUnit = "sqmeters")


# Check PA distribution from core for ablation

boxplot(HWCpoints$dist_to_core, met2_PABpoints$dist_to_core, met1_PABpoints$dist_to_core, 
        met3_PABpoints$dist_to_core, names = c("Presence", "Random PA", "3km PA", "5km PA"), 
        ylab = "Distance to core area (m)", ylab = "Dataset")

# Shows that the PA points are roughly equivalent (5km has fewer on the lower end)

# Show also in density plot
plot_df <- data.frame(Distance = c(HWCpoints$dist_to_core, met2_PABpoints$dist_to_core, met1_PABpoints$dist_to_core, 
                                   met3_PABpoints$dist_to_core),
                      Dataset = c(rep("Presence", length(HWCpoints$dist_to_core)),
                                  rep("Random PA", length(met2_PABpoints$dist_to_core)),
                                  rep("3km PA", length(met1_PABpoints$dist_to_core)),
                                  rep("5km PA", length(met3_PABpoints$dist_to_core))))
# There is suffient overlap - not complete spatial separation

# ggplots for sharing
library(ggplot2)
ggplot(plot_df, aes(Dataset, Distance)) + 
  geom_boxplot() +
  theme_bw() +
  ylab("Distance to core bear area (m)")

ggplot(plot_df,
       aes(x = Distance,
           colour = Dataset)) +
  geom_density(size = 1)

##################################
## Split into test and training ##
##################################

met1_PABpoints@data$target <- FALSE
met2_PABpoints@data$target <- FALSE
met3_PABpoints@data$target <- FALSE
HWCpoints@data$target <- TRUE

# Number of clusters
k <- 5

# Initialize fold column
HWCpoints@data$fold <- NA
set.seed(20120)
for(ds in unique(HWCpoints@data$Source)) {
  idx <- which(HWCpoints@data$Source == ds)
  HWCpoints@data$fold[idx] <- sample(
    rep(1:k, length.out = length(idx))
  )
}

table(HWCpoints@data$fold, HWCpoints@data$Source)

# Now for the PAB
set.seed(20120)
met1_PABpoints@data$fold <- sample(rep(1:k, length.out = nrow(met1_PABpoints)))
set.seed(20120)
met2_PABpoints@data$fold <- sample(rep(1:k, length.out = nrow(met2_PABpoints)))
set.seed(20120)
met3_PABpoints@data$fold <- sample(rep(1:k, length.out = nrow(met3_PABpoints)))

# Now combine the datasets for each of the cross validations
POINTS1 <- rbind(HWCpoints[c("edges_persqkm","habitat_suitability","dist_to_core", "bearloc", "target", "fold")],met1_PABpoints[c("edges_persqkm","habitat_suitability","dist_to_core","bearloc", "target", "fold")])
POINTS2 <- rbind(HWCpoints[c("edges_persqkm","habitat_suitability","dist_to_core", "bearloc", "target", "fold")],met2_PABpoints[c("edges_persqkm","habitat_suitability","dist_to_core","bearloc", "target", "fold")])
POINTS3 <- rbind(HWCpoints[c("edges_persqkm","habitat_suitability","dist_to_core", "bearloc", "target", "fold")],met3_PABpoints[c("edges_persqkm","habitat_suitability","dist_to_core","bearloc", "target", "fold")])

# For the model training, needed as factors
POINTS1@data$target <- as.factor(POINTS1@data$target)
POINTS2@data$target <- as.factor(POINTS2@data$target)
POINTS3@data$target <- as.factor(POINTS3@data$target)

#### K FOLD LOOP #####
library(pROC)

# Create a dataframe to store the results
# M1 stands for the PAB method and _1 stands for the hyperparameter methods
results <- data.frame("RF1_1" = rep(NA, 5), "RF1_2" = rep(NA, 5), "LR1" = rep(NA, 5),"GBM1" = rep(NA, 5),
                      "RF2_1" = rep(NA, 5), "RF2_2" = rep(NA, 5), "LR2" = rep(NA, 5),"GBM2" = rep(NA, 5),
                      "RF3_1" = rep(NA, 5), "RF3_2" = rep(NA, 5), "LR3" = rep(NA, 5),"GBM3" = rep(NA, 5))
rownames(results) <- c("K1", "K2", "K3", "K4", "K5")
AUC <- results
TSS <- results
AUCsd <- results
TSSsd <- results
threshold <- results

library(randomForest)
library(gbm)

all_observationsM1_2 <- c()
all_predictionsM1_2 <- c() 

# Start the loops for the analysis

for (kfold in c(1:k)){
  for (met in c(1:3)){
    spdf <- get(paste0("POINTS", met))
    # Define the training and test datasets
    # Test is fold k
    test <- spdf[spdf@data$fold == kfold,]
    training <- spdf[spdf@data$fold != kfold,]
    
    ## FIRST THE RF MODEL ##
    
    # Gen the models
    training <- training[apply(training@data,1,function(x){!any(is.na(x))}),] # remove NA's
    model1 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"])
    model2 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"], maxnodes=25)
    
    # Now predict against the test data and record key variables
    predict1 <- predict(model1, test, type = "prob")
    predict2 <- predict(model2, test, type = "prob")
    
    # For creating the pooled ROC for the end paper
    if(met == 1){
      all_observationsM1_2 <- c(all_observationsM1_2, test$target)
      all_predictionsM1_2 <- c(all_predictionsM1_2, predict2[,2])
    }
    # Record all the metrics for measuring the predictive capacity of the model
    
    # AUC (Area under ROC Curve)
    # Measures the ability of the model to discriminate between presence and absence points
    # 1 indicates perfect discrimination
    actual <- as.numeric(test@data$target)-1
    roc_curve1 <- roc(response = actual, predictor = as.numeric(predict1[,"TRUE"]))
    roc_curve2 <- roc(response = actual, predictor = as.numeric(predict2[,"TRUE"]))
    # Calculate AUC
    auc_value1 <- auc(roc_curve1)
    auc_value2 <- auc(roc_curve2)
    
    # Record the AUC
    AUC[kfold,((met*4)-3)] <- auc_value1
    AUC[kfold,((met*4)-2)] <- auc_value2
    
    # TSS (True Skill Statistic)
    # Balances sensitivity and specificity, providing a metric that is independent of prevalence
    # 1 indicates perfect agreement
    # TSS=(Sensitivity)+(Specificity)−1
    
    optimal_threshold <- coords(roc_curve1, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- optimal_threshold$sensitivity[1]
    specificity <- optimal_threshold$specificity[1]
    # Calculate TSS
    tss1 <- (sensitivity + specificity) - 1
    # Save threshold as well
    threshold[kfold,((met*4)-3)] <- optimal_threshold$threshold[1]
    
    # Same for model2
    optimal_threshold2 <- coords(roc_curve2, "best", ret = c("threshold", "sensitivity", "specificity"))
    sensitivity <- optimal_threshold2$sensitivity[1]
    specificity <- optimal_threshold2$specificity[1]
    tss2 <- (sensitivity + specificity) - 1
    # Save threshold as well
    threshold[kfold,((met*4)-2)] <- optimal_threshold2$threshold[1]
    
    # Assign to the data frame
    TSS[kfold,((met*4)-3)] <- tss1
    TSS[kfold,((met*4)-2)] <- tss2
    
    ## THEN LOGISTIC REGRESSION ##
    # First gen the model
    LR <- glm(
      target ~ .,data=training@data[, !names(training@data) %in% "fold"],
      family = binomial(link = "logit")
    )
    # Then predict based on the model
    predictions <- predict(LR,newdata=test,type="response")
    # Calculate the AUC
    roc_ob <- roc(test$target, predictions)
    aucLR <- auc(roc_ob)
    # Now TSS
    opt_threshold <- coords(roc_ob, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_threshold$sensitivity[1]
    specificity <- opt_threshold$specificity[1]
    # Calculate TSS
    tssLR <- (sensitivity + specificity) - 1
    # Save stats
    AUC[kfold,((met*4)-1)] <- aucLR
    TSS[kfold,((met*4)-1)] <- tssLR
    threshold[kfold,((met*4)-1)] <- opt_threshold$threshold[1]
    
    ## THEN GBM (GRADIENT BOOSTING)
    training@data$target <- as.numeric(training@data$target)-1
    gbm <- gbm(target ~ .,data=training@data[, !names(training@data) %in% "fold"], distribution = "bernoulli", n.trees = 500, interaction.depth = 3,
               shrinkage = 0.05, bag.fraction = 0.8)
    # Predict
    gbm_pred <- predict(gbm,newdata = test, n.trees = 500, type = "response")
    # Calculate the AUC
    roc_gbm <- roc(test$target, gbm_pred)
    aucGBM <- auc(roc_gbm)
    # Now TSS
    opt_thresh <- coords(roc_gbm, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_thresh$sensitivity[1]
    specificity <- opt_thresh$specificity[1]
    # Calculate TSS
    tssGBM <- (sensitivity + specificity) - 1
    # Save stats
    AUC[kfold,(met*4)] <- aucGBM
    TSS[kfold,(met*4)] <- tssGBM
    threshold[kfold,(met*4)] <- opt_thresh$threshold[1]
    
  }
}

################################################################################

##############################
# Spatial splitting of k-folds
###################
library(sf)

POINTS1_sf <- st_as_sf(POINTS1)
POINTS2_sf <- st_as_sf(POINTS2)
POINTS3_sf <- st_as_sf(POINTS3)

# They also need to be transformed for cv_blocks
POINTS1_sf_UTM <- st_transform(POINTS1_sf, crs = 32633)
POINTS2_sf_UTM <- st_transform(POINTS2_sf, crs = 32633)
POINTS3_sf_UTM <- st_transform(POINTS3_sf, crs = 32633)

library(blockCV) 

# Determine the value for theRange
library(spdep)
coords <- st_coordinates(POINTS1_sf_UTM)
values <- as.numeric(POINTS1_sf_UTM$target)
# Create neighbors within a distance threshold
nb <- dnearneigh(coords, 0, 70000)  # This is the end value used (0.054 Moran's I)
# and it works better with the even splitting
listw <- nb2listw(nb, style = "W")
moran.test(values, listw)

cv_blocks1 <- spatialBlock(
  speciesData = POINTS1_sf_UTM,
  species = "target",
  theRange = 70000,    # distance in meters
  selection = "random",
  k = 3
)
# Extract folds
fold_ids <- cv_blocks1$folds
# Now assign the fold IDs to the POINTS1
sp_POINTS1 <- POINTS1
for (i in seq_along(cv_blocks1$folds)) {
  test_idx <- c(cv_blocks1$folds[[i]][[2]])
  sp_POINTS1$fold[test_idx] <- i
}

#Now do for the other two
# Method 2
coords <- st_coordinates(POINTS2_sf_UTM)
values <- as.numeric(POINTS2_sf_UTM$target)
# Create neighbors within a distance threshold
nb <- dnearneigh(coords, 0, 70000)  # This is the end value used (0.054 Moran's I)
# and it works better with the even splitting
listw <- nb2listw(nb, style = "W")
moran.test(values, listw)

cv_blocks2 <- spatialBlock(
  speciesData = POINTS2_sf_UTM,
  species = "target",
  theRange = 70000,    # distance in meters
  selection = "random",
  k = 3
)
# Extract folds
fold_ids <- cv_blocks2$folds
# Now assign the fold IDs to the POINTS1
sp_POINTS2 <- POINTS2
for (i in seq_along(cv_blocks2$folds)) {
  test_idx <- c(cv_blocks2$folds[[i]][[2]])
  sp_POINTS2$fold[test_idx] <- i
}

# Method 3
coords <- st_coordinates(POINTS3_sf_UTM)
values <- as.numeric(POINTS3_sf_UTM$target)
# Create neighbors within a distance threshold
nb <- dnearneigh(coords, 0, 70000)  # This is the end value used (0.054 Moran's I)
# and it works better with the even splitting
listw <- nb2listw(nb, style = "W")
moran.test(values, listw)

cv_blocks3 <- spatialBlock(
  speciesData = POINTS3_sf_UTM,
  species = "target",
  theRange = 70000,    # distance in meters
  selection = "random",
  k = 3
)
# Extract folds
fold_ids <- cv_blocks3$folds
# Now assign the fold IDs to the POINTS1
sp_POINTS3 <- POINTS3
for (i in seq_along(cv_blocks3$folds)) {
  test_idx <- c(cv_blocks3$folds[[i]][[2]])
  sp_POINTS3$fold[test_idx] <- i
}


spresults <- data.frame("RF1_1" = rep(NA, 3), "RF1_2" = rep(NA, 3), "LR1" = rep(NA, 3),"GBM1" = rep(NA, 3),
                        "RF2_1" = rep(NA, 3), "RF2_2" = rep(NA, 3), "LR2" = rep(NA, 3),"GBM2" = rep(NA, 3),
                        "RF3_1" = rep(NA, 3), "RF3_2" = rep(NA, 3), "LR3" = rep(NA, 3),"GBM3" = rep(NA, 3))
rownames(spresults) <- c("K1", "K2", "K3")
spAUC <- spresults
spTSS <- spresults
spthreshold <- spresults

k <- 3

for (kfold in c(1:k)){
  for (met in c(1:3)){
    spdf <- get(paste0("sp_POINTS", met))
    # Define the training and test datasets
    # Test is fold k
    test <- spdf[spdf@data$fold == kfold,]
    training <- spdf[spdf@data$fold != kfold,]
    
    ## FIRST THE RF MODEL ##
    
    # Gen the models
    training <- training[apply(training@data,1,function(x){!any(is.na(x))}),] # remove NA's
    model1 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"])
    model2 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"], maxnodes=25)
    
    # Now predict against the test data and record key variables
    predict1 <- predict(model1, test, type = "prob")
    predict2 <- predict(model2, test, type = "prob")
    
    # For creating the pooled ROC for the end paper
    if(met == 1){
      all_observationsM1_2 <- c(all_observationsM1_2, test$target)
      all_predictionsM1_2 <- c(all_predictionsM1_2, predict2[,2])
    }
    # Record all the metrics for measuring the predictive capacity of the model
    
    # AUC (Area under ROC Curve)
    # Measures the ability of the model to discriminate between presence and absence points
    # 1 indicates perfect discrimination
    actual <- as.numeric(test@data$target)-1
    roc_curve1 <- roc(response = actual, predictor = as.numeric(predict1[,"TRUE"]))
    roc_curve2 <- roc(response = actual, predictor = as.numeric(predict2[,"TRUE"]))
    # Calculate AUC
    auc_value1 <- auc(roc_curve1)
    auc_value2 <- auc(roc_curve2)
    
    # Record the AUC
    spAUC[kfold,((met*4)-3)] <- auc_value1
    spAUC[kfold,((met*4)-2)] <- auc_value2
    
    # TSS (True Skill Statistic)
    # Balances sensitivity and specificity, providing a metric that is independent of prevalence
    # 1 indicates perfect agreement
    # TSS=(Sensitivity)+(Specificity)−1
    
    optimal_threshold <- coords(roc_curve1, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- optimal_threshold$sensitivity[1]
    specificity <- optimal_threshold$specificity[1]
    # Calculate TSS
    tss1 <- (sensitivity + specificity) - 1
    # Save threshold as well
    spthreshold[kfold,((met*4)-3)] <- optimal_threshold$threshold[1]
    
    # Same for model2
    optimal_threshold2 <- coords(roc_curve2, "best", ret = c("threshold", "sensitivity", "specificity"))
    sensitivity <- optimal_threshold2$sensitivity[1]
    specificity <- optimal_threshold2$specificity[1]
    tss2 <- (sensitivity + specificity) - 1
    # Save threshold as well
    spthreshold[kfold,((met*4)-2)] <- optimal_threshold2$threshold[1]
    
    # Assign to the data frame
    spTSS[kfold,((met*4)-3)] <- tss1
    spTSS[kfold,((met*4)-2)] <- tss2
    
    ## THEN LOGISTIC REGRESSION ##
    # First gen the model
    LR <- glm(
      target ~ .,data=training@data[, !names(training@data) %in% "fold"],
      family = binomial(link = "logit")
    )
    # Then predict based on the model
    predictions <- predict(LR,newdata=test,type="response")
    # Calculate the AUC
    roc_ob <- roc(test$target, predictions)
    aucLR <- auc(roc_ob)
    # Now TSS
    opt_threshold <- coords(roc_ob, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_threshold$sensitivity[1]
    specificity <- opt_threshold$specificity[1]
    # Calculate TSS
    tssLR <- (sensitivity + specificity) - 1
    # Save stats
    spAUC[kfold,((met*4)-1)] <- aucLR
    spTSS[kfold,((met*4)-1)] <- tssLR
    spthreshold[kfold,((met*4)-1)] <- opt_threshold$threshold[1]
    
    ## THEN GBM (GRADIENT BOOSTING)
    training@data$target <- as.numeric(training@data$target)-1
    gbm <- gbm(target ~ .,data=training@data[, !names(training@data) %in% "fold"], distribution = "bernoulli", n.trees = 500, interaction.depth = 3,
               shrinkage = 0.05, bag.fraction = 0.8)
    # Predict
    gbm_pred <- predict(gbm,newdata = test, n.trees = 500, type = "response")
    # Calculate the AUC
    roc_gbm <- roc(test$target, gbm_pred)
    aucGBM <- auc(roc_gbm)
    # Now TSS
    opt_thresh <- coords(roc_gbm, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_thresh$sensitivity[1]
    specificity <- opt_thresh$specificity[1]
    # Calculate TSS
    tssGBM <- (sensitivity + specificity) - 1
    # Save stats
    spAUC[kfold,(met*4)] <- aucGBM
    spTSS[kfold,(met*4)] <- tssGBM
    spthreshold[kfold,(met*4)] <- opt_thresh$threshold[1]
    
  }
}


#-----------------------
##### Final Model ######
#-----------------------
# From the results of the rest of the analysis below, the final model chosen is
# Model 1_2 from the first bear activity method

# Create the pooled ROC for the paper (Model 1_2)
rocpool <- roc(all_observationsM1_2, all_predictionsM1_2)

plot(rocpool)
auc(rocpool) #practically the same as averaged over the folds

#png("Figure4.png", units = "px", width = 3543, height = 3543, res = 300)
#plot(rocpool)
#dev.off()

total <- POINTS1[apply(POINTS1@data,1,function(x){!any(is.na(x))}),] # remove NA's

# Check the mtry is good
tuneRF(total@data[, !names(total@data) %in% "fold"], total$target, stepFactor = 3)
# Gen the model
model <- randomForest(target~.,data=total@data[, !names(total@data) %in% "fold"], maxnodes=25)

# Get the variable importance
importance(model)
varImpPlot(model)

model

#saveRDS(model, "dynamic_conflict_model.rds")

#############################################
## Generate map results ##
#############################################

cols <- colorRampPalette(c(HotspotColor))(50)

Janper <- getValues(core_area)
Janper[is.na(Janper)] <- monthper$Percentage[1]
Janper[Janper==1] <- 100-monthper$Percentage[1]

# Take a look at probabilities in some key months
JanresponseDF <- data.frame(edges_persqkm=values(edges_persqkm),
                            habitat_suitability=getValues(HSM),
                            dist_to_core=dist_to_core@data@values,
                            bearloc=Janper)

JanresponseDF <- cbind(JanresponseDF,predict(model, JanresponseDF,"prob"))
Jan1 <- classification
classification@data@values <- JanresponseDF$`TRUE`
classification[classification < 0.15] <- NA
#classification <- classification*0.5
leaflet() %>% addTiles() %>% 
  #addHeatmap(data=track_coordinates,max = 25, radius = 20, blur = 20) %>% 
  #addMarkers(data=HWCpoints, popup=HWCdetails$tags) %>% 
  addMeasure(primaryLengthUnit = "meters",primaryAreaUnit = "sqmeters") %>% 
  addRasterImage(x=classification, 
                 col=palette,
                 opacity = .7)
plot(classification, col = cols, zlim = c(0,1))

# Rerun for different months (change lines 938-939)
SepresponseDF <- JanresponseDF[, -(5:6)]
Sepper <- getValues(core_area)
# Change bracket number to generate different month
Sepper[is.na(Sepper)] <- monthper$Percentage[9]
Sepper[Sepper==1] <- 100-monthper$Percentage[9]
SepresponseDF$bearloc <- Sepper
SepresponseDF <- cbind(SepresponseDF,predict(model, SepresponseDF,"prob"))
classification@data@values <- SepresponseDF$`TRUE`
classification[classification < 0.15] <- NA
#classification <- classification*0.5
plot(classification, col = cols, zlim = c(0,1))

#writeRaster(classification, "Feb.tif", format = "GTIFF", overwrite = TRUE)







################################################################################
## CONFLICT ZONE ACTIVITY as temporal input
################################################################################
k <- 5
# Now combine the datasets for each of the cross validations
POINTS1 <- rbind(HWCpoints[c("edges_persqkm","habitat_suitability","dist_to_core", "bearact", "target", "fold")],met1_PABpoints[c("edges_persqkm","habitat_suitability","dist_to_core","bearact", "target", "fold")])
POINTS2 <- rbind(HWCpoints[c("edges_persqkm","habitat_suitability","dist_to_core", "bearact", "target", "fold")],met2_PABpoints[c("edges_persqkm","habitat_suitability","dist_to_core","bearact", "target", "fold")])
POINTS3 <- rbind(HWCpoints[c("edges_persqkm","habitat_suitability","dist_to_core", "bearact", "target", "fold")],met3_PABpoints[c("edges_persqkm","habitat_suitability","dist_to_core","bearact", "target", "fold")])

# For the model training, needed as factors
POINTS1@data$target <- as.factor(POINTS1@data$target)
POINTS2@data$target <- as.factor(POINTS2@data$target)
POINTS3@data$target <- as.factor(POINTS3@data$target)

#### K FOLD LOOP #####
library(pROC)

# Create a dataframe to store the results
# M1 stands for the PAB method and _1 stands for the hyperparameter methods
results <- data.frame("RF1_1" = rep(NA, 5), "RF1_2" = rep(NA, 5), "LR1" = rep(NA, 5),"GBM1" = rep(NA, 5),
                      "RF2_1" = rep(NA, 5), "RF2_2" = rep(NA, 5), "LR2" = rep(NA, 5),"GBM2" = rep(NA, 5),
                      "RF3_1" = rep(NA, 5), "RF3_2" = rep(NA, 5), "LR3" = rep(NA, 5),"GBM3" = rep(NA, 5))
rownames(results) <- c("K1", "K2", "K3", "K4", "K5")
AUC2 <- results
TSStwo <- results
AUCsd2 <- results
TSSsd2 <- results
threshold2 <- results

library(randomForest)
library(gbm)

all_observationsM1_2 <- c()
all_predictionsM1_2 <- c() 

# Start the loops for the analysis

for (kfold in c(1:k)){
  for (met in c(1:3)){
    spdf <- get(paste0("POINTS", met))
    # Define the training and test datasets
    # Test is fold k
    test <- spdf[spdf@data$fold == kfold,]
    training <- spdf[spdf@data$fold != kfold,]
    
    ## FIRST THE RF MODEL ##
    
    # Gen the models
    training <- training[apply(training@data,1,function(x){!any(is.na(x))}),] # remove NA's
    model1 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"])
    model2 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"], maxnodes=25)
    
    # Now predict against the test data and record key variables
    predict1 <- predict(model1, test, type = "prob")
    predict2 <- predict(model2, test, type = "prob")
    
    # For creating the pooled ROC for the end paper
    if(met == 1){
      all_observationsM1_2 <- c(all_observationsM1_2, test$target)
      all_predictionsM1_2 <- c(all_predictionsM1_2, predict2[,2])
    }
    # Record all the metrics for measuring the predictive capacity of the model
    
    # AUC (Area under ROC Curve)
    # Measures the ability of the model to discriminate between presence and absence points
    # 1 indicates perfect discrimination
    actual <- as.numeric(test@data$target)-1
    roc_curve1 <- roc(response = actual, predictor = as.numeric(predict1[,"TRUE"]))
    roc_curve2 <- roc(response = actual, predictor = as.numeric(predict2[,"TRUE"]))
    # Calculate AUC
    auc_value1 <- auc(roc_curve1)
    auc_value2 <- auc(roc_curve2)
    
    # Record the AUC
    AUC2[kfold,((met*4)-3)] <- auc_value1
    AUC2[kfold,((met*4)-2)] <- auc_value2
    
    # TSS (True Skill Statistic)
    # Balances sensitivity and specificity, providing a metric that is independent of prevalence
    # 1 indicates perfect agreement
    # TSS=(Sensitivity)+(Specificity)−1
    
    optimal_threshold <- coords(roc_curve1, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- optimal_threshold$sensitivity[1]
    specificity <- optimal_threshold$specificity[1]
    # Calculate TSS
    tss1 <- (sensitivity + specificity) - 1
    # Save threshold as well
    threshold2[kfold,((met*4)-3)] <- optimal_threshold$threshold[1]
    
    # Same for model2
    optimal_threshold2 <- coords(roc_curve2, "best", ret = c("threshold", "sensitivity", "specificity"))
    sensitivity <- optimal_threshold2$sensitivity[1]
    specificity <- optimal_threshold2$specificity[1]
    tss2 <- (sensitivity + specificity) - 1
    # Save threshold as well
    threshold2[kfold,((met*4)-2)] <- optimal_threshold2$threshold[1]
    
    # Assign to the data frame
    TSStwo[kfold,((met*4)-3)] <- tss1
    TSStwo[kfold,((met*4)-2)] <- tss2
    
    ## THEN LOGISTIC REGRESSION ##
    # First gen the model
    LR <- glm(
      target ~ .,data=training@data[, !names(training@data) %in% "fold"],
      family = binomial(link = "logit")
    )
    # Then predict based on the model
    predictions <- predict(LR,newdata=test,type="response")
    # Calculate the AUC
    roc_ob <- roc(test$target, predictions)
    aucLR <- auc(roc_ob)
    # Now TSS
    opt_threshold <- coords(roc_ob, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_threshold$sensitivity[1]
    specificity <- opt_threshold$specificity[1]
    # Calculate TSS
    tssLR <- (sensitivity + specificity) - 1
    # Save stats
    AUC2[kfold,((met*4)-1)] <- aucLR
    TSStwo[kfold,((met*4)-1)] <- tssLR
    threshold2[kfold,((met*4)-1)] <- opt_threshold$threshold[1]
    
    ## THEN GBM (GRADIENT BOOSTING)
    training@data$target <- as.numeric(training@data$target)-1
    gbm <- gbm(target ~ .,data=training@data[, !names(training@data) %in% "fold"], distribution = "bernoulli", n.trees = 500, interaction.depth = 3,
               shrinkage = 0.05, bag.fraction = 0.8)
    # Predict
    gbm_pred <- predict(gbm,newdata = test, n.trees = 500, type = "response")
    # Calculate the AUC
    roc_gbm <- roc(test$target, gbm_pred)
    aucGBM <- auc(roc_gbm)
    # Now TSS
    opt_thresh <- coords(roc_gbm, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_thresh$sensitivity[1]
    specificity <- opt_thresh$specificity[1]
    # Calculate TSS
    tssGBM <- (sensitivity + specificity) - 1
    # Save stats
    AUC2[kfold,(met*4)] <- aucGBM
    TSStwo[kfold,(met*4)] <- tssGBM
    threshold2[kfold,(met*4)] <- opt_thresh$threshold[1]
    
  }
}

##############################
# Spatial splitting of k-folds
###################
library(sf)

POINTS1_sf <- st_as_sf(POINTS1)
POINTS2_sf <- st_as_sf(POINTS2)
POINTS3_sf <- st_as_sf(POINTS3)

# They also need to be transformed for cv_blocks
POINTS1_sf_UTM <- st_transform(POINTS1_sf, crs = 32633)
POINTS2_sf_UTM <- st_transform(POINTS2_sf, crs = 32633)
POINTS3_sf_UTM <- st_transform(POINTS3_sf, crs = 32633)

library(blockCV) 

# Determine the value for theRange
library(spdep)
coords <- st_coordinates(POINTS1_sf_UTM)
values <- as.numeric(POINTS1_sf_UTM$target)
# Create neighbors within a distance threshold
nb <- dnearneigh(coords, 0, 70000)  # This is the end value used (0.054 Moran's I)
# and it works better with the even splitting
listw <- nb2listw(nb, style = "W")
moran.test(values, listw)

cv_blocks1 <- spatialBlock(
  speciesData = POINTS1_sf_UTM,
  species = "target",
  theRange = 70000,    # distance in meters
  selection = "random",
  k = 3
)
# Extract folds
fold_ids <- cv_blocks1$folds
# Now assign the fold IDs to the POINTS1
sp_POINTS1 <- POINTS1
for (i in seq_along(cv_blocks1$folds)) {
  test_idx <- c(cv_blocks1$folds[[i]][[2]])
  sp_POINTS1$fold[test_idx] <- i
}

#Now do for the other two
# Method 2
coords <- st_coordinates(POINTS2_sf_UTM)
values <- as.numeric(POINTS2_sf_UTM$target)
# Create neighbors within a distance threshold
nb <- dnearneigh(coords, 0, 70000)  # This is the end value used (0.054 Moran's I)
# and it works better with the even splitting
listw <- nb2listw(nb, style = "W")
moran.test(values, listw)

cv_blocks2 <- spatialBlock(
  speciesData = POINTS2_sf_UTM,
  species = "target",
  theRange = 70000,    # distance in meters
  selection = "random",
  k = 3
)
# Extract folds
fold_ids <- cv_blocks2$folds
# Now assign the fold IDs to the POINTS1
sp_POINTS2 <- POINTS2
for (i in seq_along(cv_blocks2$folds)) {
  test_idx <- c(cv_blocks2$folds[[i]][[2]])
  sp_POINTS2$fold[test_idx] <- i
}

# Method 3
coords <- st_coordinates(POINTS3_sf_UTM)
values <- as.numeric(POINTS3_sf_UTM$target)
# Create neighbors within a distance threshold
nb <- dnearneigh(coords, 0, 70000)  # This is the end value used (0.054 Moran's I)
# and it works better with the even splitting
listw <- nb2listw(nb, style = "W")
moran.test(values, listw)

cv_blocks3 <- spatialBlock(
  speciesData = POINTS3_sf_UTM,
  species = "target",
  theRange = 70000,    # distance in meters
  selection = "random",
  k = 3
)
# Extract folds
fold_ids <- cv_blocks3$folds
# Now assign the fold IDs to the POINTS1
sp_POINTS3 <- POINTS3
for (i in seq_along(cv_blocks3$folds)) {
  test_idx <- c(cv_blocks3$folds[[i]][[2]])
  sp_POINTS3$fold[test_idx] <- i
}


spresults <- data.frame("RF1_1" = rep(NA, 3), "RF1_2" = rep(NA, 3), "LR1" = rep(NA, 3),"GBM1" = rep(NA, 3),
                        "RF2_1" = rep(NA, 3), "RF2_2" = rep(NA, 3), "LR2" = rep(NA, 3),"GBM2" = rep(NA, 3),
                        "RF3_1" = rep(NA, 3), "RF3_2" = rep(NA, 3), "LR3" = rep(NA, 3),"GBM3" = rep(NA, 3))
rownames(spresults) <- c("K1", "K2", "K3")
spAUC2 <- spresults
spTSS2 <- spresults
spthreshold2 <- spresults

k <- 3

all_observationsM1_2 <- c()
all_predictionsM1_2 <- c() 

for (kfold in c(1:k)){
  for (met in c(1:3)){
    spdf <- get(paste0("sp_POINTS", met))
    # Define the training and test datasets
    # Test is fold k
    test <- spdf[spdf@data$fold == kfold,]
    training <- spdf[spdf@data$fold != kfold,]
    
    ## FIRST THE RF MODEL ##
    
    # Gen the models
    training <- training[apply(training@data,1,function(x){!any(is.na(x))}),] # remove NA's
    model1 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"])
    model2 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"], maxnodes=25)
    
    # Now predict against the test data and record key variables
    predict1 <- predict(model1, test, type = "prob")
    predict2 <- predict(model2, test, type = "prob")
    
    # For creating the pooled ROC for the end paper
    if(met == 1){
      all_observationsM1_2 <- c(all_observationsM1_2, test$target)
      all_predictionsM1_2 <- c(all_predictionsM1_2, predict2[,2])
    }
    # Record all the metrics for measuring the predictive capacity of the model
    
    # AUC (Area under ROC Curve)
    # Measures the ability of the model to discriminate between presence and absence points
    # 1 indicates perfect discrimination
    actual <- as.numeric(test@data$target)-1
    roc_curve1 <- roc(response = actual, predictor = as.numeric(predict1[,"TRUE"]))
    roc_curve2 <- roc(response = actual, predictor = as.numeric(predict2[,"TRUE"]))
    # Calculate AUC
    auc_value1 <- auc(roc_curve1)
    auc_value2 <- auc(roc_curve2)
    
    # Record the AUC
    spAUC2[kfold,((met*4)-3)] <- auc_value1
    spAUC2[kfold,((met*4)-2)] <- auc_value2
    
    # TSS (True Skill Statistic)
    # Balances sensitivity and specificity, providing a metric that is independent of prevalence
    # 1 indicates perfect agreement
    # TSS=(Sensitivity)+(Specificity)−1
    
    optimal_threshold <- coords(roc_curve1, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- optimal_threshold$sensitivity[1]
    specificity <- optimal_threshold$specificity[1]
    # Calculate TSS
    tss1 <- (sensitivity + specificity) - 1
    # Save threshold as well
    spthreshold2[kfold,((met*4)-3)] <- optimal_threshold$threshold[1]
    
    # Same for model2
    optimal_threshold2 <- coords(roc_curve2, "best", ret = c("threshold", "sensitivity", "specificity"))
    sensitivity <- optimal_threshold2$sensitivity[1]
    specificity <- optimal_threshold2$specificity[1]
    tss2 <- (sensitivity + specificity) - 1
    # Save threshold as well
    spthreshold2[kfold,((met*4)-2)] <- optimal_threshold2$threshold[1]
    
    # Assign to the data frame
    spTSS2[kfold,((met*4)-3)] <- tss1
    spTSS2[kfold,((met*4)-2)] <- tss2
    
    ## THEN LOGISTIC REGRESSION ##
    # First gen the model
    LR <- glm(
      target ~ .,data=training@data[, !names(training@data) %in% "fold"],
      family = binomial(link = "logit")
    )
    # Then predict based on the model
    predictions <- predict(LR,newdata=test,type="response")
    # Calculate the AUC
    roc_ob <- roc(test$target, predictions)
    aucLR <- auc(roc_ob)
    # Now TSS
    opt_threshold <- coords(roc_ob, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_threshold$sensitivity[1]
    specificity <- opt_threshold$specificity[1]
    # Calculate TSS
    tssLR <- (sensitivity + specificity) - 1
    # Save stats
    spAUC2[kfold,((met*4)-1)] <- aucLR
    spTSS2[kfold,((met*4)-1)] <- tssLR
    spthreshold2[kfold,((met*4)-1)] <- opt_threshold$threshold[1]
    
    ## THEN GBM (GRADIENT BOOSTING)
    training@data$target <- as.numeric(training@data$target)-1
    gbm <- gbm(target ~ .,data=training@data[, !names(training@data) %in% "fold"], distribution = "bernoulli", n.trees = 500, interaction.depth = 3,
               shrinkage = 0.05, bag.fraction = 0.8)
    # Predict
    gbm_pred <- predict(gbm,newdata = test, n.trees = 500, type = "response")
    # Calculate the AUC
    roc_gbm <- roc(test$target, gbm_pred)
    aucGBM <- auc(roc_gbm)
    # Now TSS
    opt_thresh <- coords(roc_gbm, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_thresh$sensitivity[1]
    specificity <- opt_thresh$specificity[1]
    # Calculate TSS
    tssGBM <- (sensitivity + specificity) - 1
    # Save stats
    spAUC2[kfold,(met*4)] <- aucGBM
    spTSS2[kfold,(met*4)] <- tssGBM
    spthreshold2[kfold,(met*4)] <- opt_thresh$threshold[1]
    
  }
}





################################################################################
#### DISTANCE FROM CORE BEAR AREA REMOVED ######################################
################################################################################
# As part of the ablation study, we also tested the model results with the core
# bear area variable removed

k <- 5

# Now combine the datasets for each of the cross validations
POINTS1 <- rbind(HWCpoints[c("edges_persqkm","habitat_suitability","bearloc", "target", "fold")],met1_PABpoints[c("edges_persqkm","habitat_suitability","bearloc", "target", "fold")])
POINTS2 <- rbind(HWCpoints[c("edges_persqkm","habitat_suitability", "bearloc", "target", "fold")],met2_PABpoints[c("edges_persqkm","habitat_suitability","bearloc", "target", "fold")])
POINTS3 <- rbind(HWCpoints[c("edges_persqkm","habitat_suitability", "bearloc", "target", "fold")],met3_PABpoints[c("edges_persqkm","habitat_suitability","bearloc", "target", "fold")])

# For the model training, needed as factors
POINTS1@data$target <- as.factor(POINTS1@data$target)
POINTS2@data$target <- as.factor(POINTS2@data$target)
POINTS3@data$target <- as.factor(POINTS3@data$target)

#### K FOLD LOOP #####
library(pROC)

# Create a dataframe to store the results
# M1 stands for the PAB method and _1 stands for the hyperparameter methods
results <- data.frame("RF1_1" = rep(NA, 5), "RF1_2" = rep(NA, 5), "LR1" = rep(NA, 5),"GBM1" = rep(NA, 5),
                      "RF2_1" = rep(NA, 5), "RF2_2" = rep(NA, 5), "LR2" = rep(NA, 5),"GBM2" = rep(NA, 5),
                      "RF3_1" = rep(NA, 5), "RF3_2" = rep(NA, 5), "LR3" = rep(NA, 5),"GBM3" = rep(NA, 5))
rownames(results) <- c("K1", "K2", "K3", "K4", "K5")
AUC3 <- results
TSS3 <- results
AUCsd3 <- results
TSSsd3 <- results
threshold3 <- results

library(randomForest)
library(gbm)

all_observationsM1_2 <- c()
all_predictionsM1_2 <- c() 

# Start the loops for the analysis
for (kfold in c(1:k)){
  for (met in c(1:3)){
    spdf <- get(paste0("POINTS", met))
    # Define the training and test datasets
    # Test is fold k
    test <- spdf[spdf@data$fold == kfold,]
    training <- spdf[spdf@data$fold != kfold,]
    
    ## FIRST THE RF MODEL ##
    
    # Gen the models
    training <- training[apply(training@data,1,function(x){!any(is.na(x))}),] # remove NA's
    model1 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"])
    model2 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"], maxnodes=25)
    
    # Now predict against the test data and record key variables
    predict1 <- predict(model1, test, type = "prob")
    predict2 <- predict(model2, test, type = "prob")
    
    # For creating the pooled ROC for the end paper
    if(met == 1){
      all_observationsM1_2 <- c(all_observationsM1_2, test$target)
      all_predictionsM1_2 <- c(all_predictionsM1_2, predict2[,2])
    }
    # Record all the metrics for measuring the predictive capacity of the model
    
    # AUC (Area under ROC Curve)
    # Measures the ability of the model to discriminate between presence and absence points
    # 1 indicates perfect discrimination
    actual <- as.numeric(test@data$target)-1
    roc_curve1 <- roc(response = actual, predictor = as.numeric(predict1[,"TRUE"]))
    roc_curve2 <- roc(response = actual, predictor = as.numeric(predict2[,"TRUE"]))
    # Calculate AUC
    auc_value1 <- auc(roc_curve1)
    auc_value2 <- auc(roc_curve2)
    
    # Record the AUC
    AUC3[kfold,((met*4)-3)] <- auc_value1
    AUC3[kfold,((met*4)-2)] <- auc_value2
    
    # TSS (True Skill Statistic)
    # Balances sensitivity and specificity, providing a metric that is independent of prevalence
    # 1 indicates perfect agreement
    # TSS=(Sensitivity)+(Specificity)−1
    
    optimal_threshold <- coords(roc_curve1, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- optimal_threshold$sensitivity[1]
    specificity <- optimal_threshold$specificity[1]
    # Calculate TSS
    tss1 <- (sensitivity + specificity) - 1
    # Save threshold as well
    threshold3[kfold,((met*4)-3)] <- optimal_threshold$threshold[1]
    
    # Same for model2
    optimal_threshold2 <- coords(roc_curve2, "best", ret = c("threshold", "sensitivity", "specificity"))
    sensitivity <- optimal_threshold2$sensitivity[1]
    specificity <- optimal_threshold2$specificity[1]
    tss2 <- (sensitivity + specificity) - 1
    # Save threshold as well
    threshold3[kfold,((met*4)-2)] <- optimal_threshold2$threshold[1]
    
    # Assign to the data frame
    TSS3[kfold,((met*4)-3)] <- tss1
    TSS3[kfold,((met*4)-2)] <- tss2
    
    ## THEN LOGISTIC REGRESSION ##
    # First gen the model
    LR <- glm(
      target ~ .,data=training@data[, !names(training@data) %in% "fold"],
      family = binomial(link = "logit")
    )
    # Then predict based on the model
    predictions <- predict(LR,newdata=test,type="response")
    # Calculate the AUC
    roc_ob <- roc(test$target, predictions)
    aucLR <- auc(roc_ob)
    # Now TSS
    opt_threshold <- coords(roc_ob, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_threshold$sensitivity[1]
    specificity <- opt_threshold$specificity[1]
    # Calculate TSS
    tssLR <- (sensitivity + specificity) - 1
    # Save stats
    AUC3[kfold,((met*4)-1)] <- aucLR
    TSS3[kfold,((met*4)-1)] <- tssLR
    threshold3[kfold,((met*4)-1)] <- opt_threshold$threshold[1]
    
    ## THEN GBM (GRADIENT BOOSTING)
    training@data$target <- as.numeric(training@data$target)-1
    gbm <- gbm(target ~ .,data=training@data[, !names(training@data) %in% "fold"], distribution = "bernoulli", n.trees = 500, interaction.depth = 3,
               shrinkage = 0.05, bag.fraction = 0.8)
    # Predict
    gbm_pred <- predict(gbm,newdata = test, n.trees = 500, type = "response")
    # Calculate the AUC
    roc_gbm <- roc(test$target, gbm_pred)
    aucGBM <- auc(roc_gbm)
    # Now TSS
    opt_thresh <- coords(roc_gbm, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_thresh$sensitivity[1]
    specificity <- opt_thresh$specificity[1]
    # Calculate TSS
    tssGBM <- (sensitivity + specificity) - 1
    # Save stats
    AUC3[kfold,(met*4)] <- aucGBM
    TSS3[kfold,(met*4)] <- tssGBM
    threshold3[kfold,(met*4)] <- opt_thresh$threshold[1]
    
  }
}


##############################
# Spatial splitting of k-folds
###################
library(sf)

POINTS1_sf <- st_as_sf(POINTS1)
POINTS2_sf <- st_as_sf(POINTS2)
POINTS3_sf <- st_as_sf(POINTS3)

# They also need to be transformed for cv_blocks
POINTS1_sf_UTM <- st_transform(POINTS1_sf, crs = 32633)
POINTS2_sf_UTM <- st_transform(POINTS2_sf, crs = 32633)
POINTS3_sf_UTM <- st_transform(POINTS3_sf, crs = 32633)

library(blockCV) 

# Determine the value for theRange
library(spdep)
coords <- st_coordinates(POINTS1_sf_UTM)
values <- as.numeric(POINTS1_sf_UTM$target)
# Create neighbors within a distance threshold
nb <- dnearneigh(coords, 0, 70000)  # This is the end value used (0.054 Moran's I)
# and it works better with the even splitting
listw <- nb2listw(nb, style = "W")
moran.test(values, listw)

cv_blocks1 <- spatialBlock(
  speciesData = POINTS1_sf_UTM,
  species = "target",
  theRange = 70000,    # distance in meters
  selection = "random",
  k = 3
)
# Extract folds
fold_ids <- cv_blocks1$folds
# Now assign the fold IDs to the POINTS1
sp_POINTS1 <- POINTS1
for (i in seq_along(cv_blocks1$folds)) {
  test_idx <- c(cv_blocks1$folds[[i]][[2]])
  sp_POINTS1$fold[test_idx] <- i
}

#Now do for the other two
# Method 2
coords <- st_coordinates(POINTS2_sf_UTM)
values <- as.numeric(POINTS2_sf_UTM$target)
# Create neighbors within a distance threshold
nb <- dnearneigh(coords, 0, 70000)  # This is the end value used (0.054 Moran's I)
# and it works better with the even splitting
listw <- nb2listw(nb, style = "W")
moran.test(values, listw)

cv_blocks2 <- spatialBlock(
  speciesData = POINTS2_sf_UTM,
  species = "target",
  theRange = 70000,    # distance in meters
  selection = "random",
  k = 3
)
# Extract folds
fold_ids <- cv_blocks2$folds
# Now assign the fold IDs to the POINTS1
sp_POINTS2 <- POINTS2
for (i in seq_along(cv_blocks2$folds)) {
  test_idx <- c(cv_blocks2$folds[[i]][[2]])
  sp_POINTS2$fold[test_idx] <- i
}

# Method 3
coords <- st_coordinates(POINTS3_sf_UTM)
values <- as.numeric(POINTS3_sf_UTM$target)
# Create neighbors within a distance threshold
nb <- dnearneigh(coords, 0, 70000)  # This is the end value used (0.054 Moran's I)
# and it works better with the even splitting
listw <- nb2listw(nb, style = "W")
moran.test(values, listw)

cv_blocks3 <- spatialBlock(
  speciesData = POINTS3_sf_UTM,
  species = "target",
  theRange = 70000,    # distance in meters
  selection = "random",
  k = 3
)
# Extract folds
fold_ids <- cv_blocks3$folds
# Now assign the fold IDs to the POINTS1
sp_POINTS3 <- POINTS3
for (i in seq_along(cv_blocks3$folds)) {
  test_idx <- c(cv_blocks3$folds[[i]][[2]])
  sp_POINTS3$fold[test_idx] <- i
}


spresults <- data.frame("RF1_1" = rep(NA, 3), "RF1_2" = rep(NA, 3), "LR1" = rep(NA, 3),"GBM1" = rep(NA, 3),
                        "RF2_1" = rep(NA, 3), "RF2_2" = rep(NA, 3), "LR2" = rep(NA, 3),"GBM2" = rep(NA, 3),
                        "RF3_1" = rep(NA, 3), "RF3_2" = rep(NA, 3), "LR3" = rep(NA, 3),"GBM3" = rep(NA, 3))
rownames(spresults) <- c("K1", "K2", "K3")
spAUC3 <- spresults
spTSS3 <- spresults
spthreshold3 <- spresults

k <- 3

for (kfold in c(1:k)){
  for (met in c(1:3)){
    spdf <- get(paste0("sp_POINTS", met))
    # Define the training and test datasets
    # Test is fold k
    test <- spdf[spdf@data$fold == kfold,]
    training <- spdf[spdf@data$fold != kfold,]
    
    ## FIRST THE RF MODEL ##
    
    # Gen the models
    training <- training[apply(training@data,1,function(x){!any(is.na(x))}),] # remove NA's
    model1 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"])
    model2 <- randomForest(target~.,data=training@data[, !names(training@data) %in% "fold"], maxnodes=25)
    
    # Now predict against the test data and record key variables
    predict1 <- predict(model1, test, type = "prob")
    predict2 <- predict(model2, test, type = "prob")
    
    # For creating the pooled ROC for the end paper
    if(met == 1){
      all_observationsM1_2 <- c(all_observationsM1_2, test$target)
      all_predictionsM1_2 <- c(all_predictionsM1_2, predict2[,2])
    }
    # Record all the metrics for measuring the predictive capacity of the model
    
    # AUC (Area under ROC Curve)
    # Measures the ability of the model to discriminate between presence and absence points
    # 1 indicates perfect discrimination
    actual <- as.numeric(test@data$target)-1
    roc_curve1 <- roc(response = actual, predictor = as.numeric(predict1[,"TRUE"]))
    roc_curve2 <- roc(response = actual, predictor = as.numeric(predict2[,"TRUE"]))
    # Calculate AUC
    auc_value1 <- auc(roc_curve1)
    auc_value2 <- auc(roc_curve2)
    
    # Record the AUC
    spAUC3[kfold,((met*4)-3)] <- auc_value1
    spAUC3[kfold,((met*4)-2)] <- auc_value2
    
    # TSS (True Skill Statistic)
    # Balances sensitivity and specificity, providing a metric that is independent of prevalence
    # 1 indicates perfect agreement
    # TSS=(Sensitivity)+(Specificity)−1
    
    optimal_threshold <- coords(roc_curve1, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- optimal_threshold$sensitivity[1]
    specificity <- optimal_threshold$specificity[1]
    # Calculate TSS
    tss1 <- (sensitivity + specificity) - 1
    # Save threshold as well
    spthreshold3[kfold,((met*4)-3)] <- optimal_threshold$threshold[1]
    
    # Same for model2
    optimal_threshold2 <- coords(roc_curve2, "best", ret = c("threshold", "sensitivity", "specificity"))
    sensitivity <- optimal_threshold2$sensitivity[1]
    specificity <- optimal_threshold2$specificity[1]
    tss2 <- (sensitivity + specificity) - 1
    # Save threshold as well
    spthreshold3[kfold,((met*4)-2)] <- optimal_threshold2$threshold[1]
    
    # Assign to the data frame
    spTSS3[kfold,((met*4)-3)] <- tss1
    spTSS3[kfold,((met*4)-2)] <- tss2
    
    ## THEN LOGISTIC REGRESSION ##
    # First gen the model
    LR <- glm(
      target ~ .,data=training@data[, !names(training@data) %in% "fold"],
      family = binomial(link = "logit")
    )
    # Then predict based on the model
    predictions <- predict(LR,newdata=test,type="response")
    # Calculate the AUC
    roc_ob <- roc(test$target, predictions)
    aucLR <- auc(roc_ob)
    # Now TSS
    opt_threshold <- coords(roc_ob, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_threshold$sensitivity[1]
    specificity <- opt_threshold$specificity[1]
    # Calculate TSS
    tssLR <- (sensitivity + specificity) - 1
    # Save stats
    spAUC3[kfold,((met*4)-1)] <- aucLR
    spTSS3[kfold,((met*4)-1)] <- tssLR
    spthreshold3[kfold,((met*4)-1)] <- opt_threshold$threshold[1]
    
    ## THEN GBM (GRADIENT BOOSTING)
    training@data$target <- as.numeric(training@data$target)-1
    gbm <- gbm(target ~ .,data=training@data[, !names(training@data) %in% "fold"], distribution = "bernoulli", n.trees = 500, interaction.depth = 3,
               shrinkage = 0.05, bag.fraction = 0.8)
    # Predict
    gbm_pred <- predict(gbm,newdata = test, n.trees = 500, type = "response")
    # Calculate the AUC
    roc_gbm <- roc(test$target, gbm_pred)
    aucGBM <- auc(roc_gbm)
    # Now TSS
    opt_thresh <- coords(roc_gbm, "best", ret = c("threshold", "sensitivity", "specificity"))
    # Extract sensitivity and specificity at the optimal threshold
    sensitivity <- opt_thresh$sensitivity[1]
    specificity <- opt_thresh$specificity[1]
    # Calculate TSS
    tssGBM <- (sensitivity + specificity) - 1
    # Save stats
    spAUC3[kfold,(met*4)] <- aucGBM
    spTSS3[kfold,(met*4)] <- tssGBM
    spthreshold3[kfold,(met*4)] <- opt_thresh$threshold[1]
    
  }
}


################################################################################
# Export the gathered AUC and TSS data

library(openxlsx)

# Just make a big workbook to export everything and then work there
wb <- createWorkbook()

addWorksheet(wb, "AUC")
writeData(wb, "AUC", AUC)

addWorksheet(wb, "AUC2")
writeData(wb, "AUC2", AUC2)

addWorksheet(wb, "AUC3")
writeData(wb, "AUC3", AUC3)

addWorksheet(wb, "spAUC")
writeData(wb, "spAUC", spAUC)

addWorksheet(wb, "spAUC2")
writeData(wb, "spAUC2", spAUC2)

addWorksheet(wb, "spAUC3")
writeData(wb, "spAUC3", spAUC3)

addWorksheet(wb, "TSS")
writeData(wb, "TSS", TSS)

addWorksheet(wb, "TSS2")
writeData(wb, "TSS2", TSStwo)

addWorksheet(wb, "TSS3")
writeData(wb, "TSS3", TSS3)

addWorksheet(wb, "spTSS")
writeData(wb, "spTSS", spTSS)

addWorksheet(wb, "spTSS2")
writeData(wb, "spTSS2", spTSS2)

addWorksheet(wb, "spTSS3")
writeData(wb, "spTSS3", spTSS3)

addWorksheet(wb, "threshold1")
writeData(wb, "threshold1", threshold)

addWorksheet(wb, "threshold2")
writeData(wb, "threshold2", threshold2)

addWorksheet(wb, "threshold3")
writeData(wb, "threshold3", threshold3)

addWorksheet(wb, "spthreshold")
writeData(wb, "spthreshold", spthreshold)

addWorksheet(wb, "spthreshold2")
writeData(wb, "spthreshold2", spthreshold2)

addWorksheet(wb, "spthreshold3")
writeData(wb, "spthreshold3", spthreshold3)

saveWorkbook(wb,
             "Review_Model_Results.xlsx",
             overwrite = TRUE)




