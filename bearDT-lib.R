calcBearea <- function(BearLoc,
                       time_now=Sys.time(), # current time
                       forecast=48, # number of hours for which the risk is forcasted
                       Altitude_map,
                       mean_dispersal_coef_a = 1.08046,   # intersect of dispersal model for mean dispersal
                       mean_dispersal_coef_b = 0.48566,   # coefficient of dispersal model for mean dispersal
                       upper_dispersal_coef_a = 2.0892,  # intersect of dispersal model for maximum dispersal
                       upper_dispersal_coef_b = 0.4872   # coefficient of dispersal model for maximum dispersal
) {
  BearAlt <- extract(Altitude_map,BearLoc)
  
  dtime <- as.numeric(difftime(time_now,BearLoc$datetime,units="hours"))
  dtime <- dtime + forecast
  b <- 1/0.04474
  a <- 0.71022*b
  alpha <- 2.0892
  beta <- 0.4872
  dalt <- exp(5.15398+0.30816*log(dtime))
  distance <- a*dtime/(b+dtime)
  distance[dtime<24] <- exp(alpha+beta*log(dtime/24))[dtime<24] # use exponential dispersal model in first 24 hours
  distance <- distance*1000
  distance_RASTER <- distanceFromPoints(Altitude_map,BearLoc[1,])
  distance_MASK <- distance_RASTER<distance[1]
  altitude_MASK <- Altitude_map>(BearAlt[1]-dalt[1])&Altitude_map<(BearAlt[1]+dalt[1])
  dist_alt_MASK <- altitude_MASK&distance_MASK
  #dist_alt_MASK[dist_alt_MASK==0] <- NA
  
  if(nrow(BearLoc)>1) {
    for (i in 2:nrow(BearLoc)) {
      distance_RASTER <- distanceFromPoints(Altitude_map,BearLoc[i,])
      dist <- ifelse(length(dtime)>1,distance[i],distance)
      distance_MASK <- distance_RASTER<dist
      da <- ifelse(length(dtime)>1,dalt[i],dalt)
      altitude_MASK <- Altitude_map>(BearAlt[i]-da)&Altitude_map<(BearAlt[i]+da)
      MASK <- altitude_MASK&distance_MASK
      dist_alt_MASK <- dist_alt_MASK | MASK
    }
  }
  
  # control for capes - calc distance within raster when bear has to go around a mountain or valley
  prob_dist_map <- dist_alt_MASK ==2
  for (i in 1:nrow(BearLoc)) {
    index <- cellFromXY(dist_alt_MASK, st_coordinates(BearLoc)[i,])
    dist_alt_MASK[index] <- 2
    dist_in_altrange <- try(raster::gridDistance(dist_alt_MASK,2,omit=0),silent=TRUE)
    mu <- mean_dispersal_coef_a+mean_dispersal_coef_b*log(dtime[i]/24)
    sd <- (upper_dispersal_coef_a+upper_dispersal_coef_b*log(dtime[i]/24)-mu)/2
    
    if (class(dist_in_altrange) == "try-error") { # distance covers only 1 cell
      prob_map <- dist_alt_MASK == 2
      
    } else { # distance covers multiple cells
      
      prob_map <- dist_in_altrange
      prob_map@data@values <- pnorm(log(prob_map@data@values/1000), mu, sd, lower.tail = F)
      prob_map[is.na(prob_map),] <- 0
      
    }
    dist_alt_MASK[index] <- 1 # reset value
    
    prob_dist_map <- 1-(1-prob_dist_map)*(1-prob_map)
  }
  return(prob_dist_map)
}

getHWCobservations <- function(cookie,
                               from = Sys.Date()-14) { # by default: look two weeks back in time
  # human-bear-conflict
  df <- sensingcluesr::get_observations(cookie, 
                                        group = c(), # check any group that the user has access to
                                        from = from,
                                        bounds = list(north = 43.5, east = 26, south = 42, west = 23.5), # filter on Stara Planina
                                        allAttributes = TRUE,
                                        filteredConcepts = c("https://sensingclues.poolparty.biz/SCCSSOntology/106")) #hwc
  
  # filter out the hwc observations related to bears
  HWC <- df %>% dplyr::filter(grepl("bear", tags)) %>% 
    dplyr::filter(!duplicated(entityId)) %>% 
    dplyr::mutate(datetime = lubridate::as_datetime(when, tz = Sys.timezone()), # format="%Y-%m-%dT%H:%M:%OS"
                  lon = longitude,
                  lat = latitude) %>% 
    dplyr::arrange(datetime) %>% 
    dplyr::select(any_of(c("lon","lat","entityId","datetime","longitude","latitude","tags")))
  
  # make the HWC dataframe spatial
  HWC <- st_as_sf(HWC, coords = c("lon", "lat"), crs = 4326)
  
  return(HWC)
}

calcWhodunnit <- function(suspects,
                         Bear_names,
                         incident,
                         Altitude_map,
                         mean_dispersal_coef_a = 1.08046,   # intersect of dispersal model for mean dispersal
                         mean_dispersal_coef_b = 0.48566,   # coefficient of dispersal model for mean dispersal
                         upper_dispersal_coef_a = 2.0892,  # intersect of dispersal model for maximum dispersal
                         upper_dispersal_coef_b = 0.4872   # coefficient of dispersal model for maximum dispersal
) {
  t_inc2 <- incident$datetime
  # calc probability that incident is caused by any known bears:
  dtime <- as.numeric(difftime(t_inc2,suspects$datetime,units="hours"))
  dtime[dtime<0] <- 0 # avoid math errors - also, suspects can only be calculated from incidents prior to event
  if(any(dtime<1&dtime>0)) {message("Warning: cannot calculate range for dtime less than 1 hour. Setting dtime to 1 hour.")}
  dtime[dtime<1&dtime>0] <- 1
  
  BearAlt <- extract(Altitude_map,suspects)
  
  b <- 1/0.04474
  a <- 0.71022*b
  alpha <- 2.0892
  beta <- 0.4872
  dalt <- exp(5.15398+0.30816*log(dtime))
  distance <- a*dtime/(b+dtime)
  distance[dtime<24] <- exp(alpha+beta*log(dtime/24))[dtime<24] # use exponential dispersal model in first 24 hours
  distance <- distance*1000
  distance_RASTER <- distanceFromPoints(Altitude_map,suspects[1,])
  distance_MASK <- distance_RASTER<distance[1]
  altitude_MASK <- Altitude_map>(BearAlt[1]-dalt[1])&Altitude_map<(BearAlt[1]+dalt[1])
  dist_alt_MASK <- altitude_MASK&distance_MASK
  #dist_alt_MASK[dist_alt_MASK==0] <- NA
  
  if(nrow(suspects)>1) {
    for (i in 2:nrow(suspects)) {
      distance_RASTER <- distanceFromPoints(Altitude_map,suspects[i,])
      dist <- ifelse(length(dtime)>1,distance[i],distance)
      distance_MASK <- distance_RASTER<dist
      da <- ifelse(length(dtime)>1,dalt[i],dalt)
      altitude_MASK <- Altitude_map>(BearAlt[i]-da)&Altitude_map<(BearAlt[i]+da)
      MASK <- altitude_MASK&distance_MASK
      dist_alt_MASK <- dist_alt_MASK | MASK
    }
  }
  
  # control for capes - calc distance within raster when bear has to go around a mountain or valley
  p.distance <- NULL
  prob_dist_map <- dist_alt_MASK ==2
  for (i in 1:nrow(suspects)) {
    index <- cellFromXY(dist_alt_MASK, st_coordinates(suspects)[i,])
    dist_alt_MASK[index] <- 2
    dist_in_altrange <- try(raster::gridDistance(dist_alt_MASK,2,omit=0),silent=TRUE)
    mu <- mean_dispersal_coef_a+mean_dispersal_coef_b*log(dtime[i]/24)
    sd <- (upper_dispersal_coef_a+upper_dispersal_coef_b*log(dtime[i]/24)-mu)/2
    
    if (class(dist_in_altrange) == "try-error") { # distance covers only 1 cell
      prob_map <- dist_alt_MASK == 2
      
    } else { # distance covers multiple cells
      
      prob_map <- dist_in_altrange
      prob_map@data@values <- pnorm(log(prob_map@data@values/1000), mu, sd, lower.tail = F)
      prob_map[is.na(prob_map),] <- 0
      
    }
    dist_alt_MASK[index] <- 1 # reset value
    
    pval <- extract(prob_map,incident)
    prob_dist_map <- 1-(1-prob_dist_map)*(1-prob_map)
    
    p.distance <- c(p.distance,pval)
  }
  
  p <- 1-prod(1-p.distance)
  new.bear <- p<0.05 # if p is lower than 0.05, then there is a new bear
  if(new.bear) {
    troublemaker <- "unknown"
  } else {
    # which bear was the most likely problem bear for the new incident?
    replace <- which(p.distance==max(p.distance))
    troublemaker <- Bear_names[replace]
  }
  
  return(troublemaker)
}

calcActivityArea <- function(BearLoc,dtime=1,BearAlt,bounds="median",Altitude_map) {
  if(any(dtime<1)) {message("Warning: cannot calculate range for dtime less than 1 hour. Setting dtime to 1 hour.")}
  dtime[dtime<1] <- 1
  if(bounds=="upper"){
    b <- 1/0.04474
    a <- 0.71022*b
    alpha <- 2.0892
    beta <- 0.4872
    dalt <- exp(5.15398+0.30816*log(dtime))
  } else {
    b <- 1/0.05528
    a <- 0.30425*b
    alpha <- 1.08046
    beta <- 0.48566
    dalt <- exp(3.91235+0.40901*log(dtime))
  }
  distance <- a*dtime/(b+dtime)
  distance[dtime<24] <- exp(alpha+beta*log(dtime/24))[dtime<24] # use exponential dispersal model in first 24 hours
  distance <- distance*1000
  distance_RASTER <- distanceFromPoints(Altitude_map,BearLoc[1,])
  distance_MASK <- distance_RASTER<distance[1]
  altitude_MASK <- Altitude_map>(BearAlt[1]-dalt[1])&Altitude_map<(BearAlt[1]+dalt[1])
  dist_alt_MASK <- altitude_MASK&distance_MASK
  dist_alt_MASK[dist_alt_MASK==0] <- NA
  
  if(nrow(BearLoc)>1) {
    for (i in 2:nrow(BearLoc)) {
      distance_RASTER <- distanceFromPoints(Altitude_map,BearLoc[i,])
      dist <- ifelse(length(dtime)>1,distance[i],distance)
      distance_MASK <- distance_RASTER<dist
      da <- ifelse(length(dtime)>1,dalt[i],dalt)
      altitude_MASK <- Altitude_map>(BearAlt[i]-da)&Altitude_map<(BearAlt[i]+da)
      MASK <- altitude_MASK&distance_MASK
      dist_alt_MASK <- dist_alt_MASK | MASK
    }
  }
  
  # control for capes - calc distance within raster when bear has to go around a mountain or valley
  index <- cellFromXY(dist_alt_MASK, st_coordinates(BearLoc)[1,])
  dist_alt_MASK[index] <- 2
  dist_in_altrange <- try(raster::gridDistance(dist_alt_MASK,2,omit=NA),silent=TRUE)
  if (class(dist_in_altrange) == "try-error") {
    dist_in_altrange2 <- dist_alt_MASK == 2
  } else {
    dist_in_altrange2 <- dist_in_altrange<distance[1]
  }
  
  if(length(dtime)>1) {
    for (i in 2:length(dtime)) {
      index <- cellFromXY(dist_alt_MASK, st_coordinates(BearLoc)[i,])
      dist_alt_MASK[index] <- 2
      dist_in_altrange <- try(raster::gridDistance(dist_alt_MASK,2,omit=NA),silent=TRUE)
      if (class(dist_in_altrange) == "try-error") {
        MASK <- dist_alt_MASK == 2
      } else {
        MASK <- dist_in_altrange<distance[i]
      }
      dist_in_altrange2 <- dist_in_altrange2 | MASK
    }
  } 
  altrange_poly <- rasterToPolygons(dist_in_altrange2,function(x){x==1},dissolve=TRUE)
  return(altrange_poly)
}