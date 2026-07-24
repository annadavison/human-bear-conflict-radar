############### Creation of the bear dispersal model ###########################
#-------------------------------------------------------------------------------
# Load in the bear telemetry files. These are excel spreadsheets containing the
# latitudes, longitudes, date/time and sometimes altitude of bear within the
# Stara Planina region. Full details on this data are within the associated
# publication
GPS_files <- list.files(pattern = ".xlsx")
GPS_files <- GPS_files[-1]
Bear_names <- strsplit(GPS_files[1],"_")[[1]][2]
for (i in 2:length(GPS_files)) {
  name <- strsplit(GPS_files[i],"_")[[1]][2]
  Bear_names <- c(Bear_names,name)
}

# process bear data:
process_beardata <- function(filename) {
  BearDATA <- read_xlsx(filename)
  colnames <- names(BearDATA)
  
  LON_index <- grep("lon",colnames,ignore.case = TRUE)
  LAT_index <- grep("lat",colnames,ignore.case = TRUE)
  UTC_index <- grep("utc",colnames,ignore.case = TRUE)
  Time_index <- grep("time",colnames,ignore.case = TRUE)
  UTC_time_index <- UTC_index[UTC_index %in% Time_index]
  ID_index <- grep("id",colnames,ignore.case = TRUE)
  Date_index <- grep("date",colnames,ignore.case = TRUE)
  UTC_date_index <- UTC_index[UTC_index %in% Date_index]
  ALT_index <- grep("alt",colnames,ignore.case = TRUE)
  if(length(ALT_index)==0){ALT_index <- grep("height",colnames,ignore.case = TRUE)}
  
  BearGPS <- data.frame(N=1:nrow(BearDATA))
  if(length(UTC_date_index)>0) {
    if(any(colnames == "UTC_Date")) {date <- BearDATA$UTC_Date} else {
      if(any(colnames == "UTC_DATE")) {date <- BearDATA$UTC_DATE} else {
        if(any(colnames == "UTC-DATE")) {date <- BearDATA$`UTC-DATE`} else {
          date <- unlist(as.vector(BearDATA[colnames[UTC_date_index]]),use.names = FALSE)
        }
      }
    }
  } else {
    date <- BearDATA$Date
  }
  date <- as.Date(date,tryFormats = c("%Y-%m-%d",
                                      "%d-%m-%Y",
                                      "%d/%m/%Y",
                                      "%m/%d/%Y",
                                      "%Y/%m/%d",
                                      "%d.%m.%Y"))
  if(length(UTC_time_index)>0) {
    if(any(colnames == "UTC-Time")){time <- BearDATA$`UTC-Time`} else {
      if(any(colnames == "UTC_TIME")){time <- BearDATA$UTC_TIME} else {
        time <- unlist(as.vector(BearDATA[colnames[UTC_time_index]]),use.names = FALSE)
      }
    }
  } else {
    if(any(colnames == "Time")){time <- BearDATA$Time}
    if(any(colnames == "Time_LTM")){time <- BearDATA$Time_LTM}
    if(any(colnames == "DateTime")){time <- BearDATA$DateTime}  
  }
  if(nchar(time[1])>8) {time <- substr(as.POSIXlt(time),12,19)} else {
    if(nchar(time[1])==5) {time <- paste(time,"00",sep=":")}
  }
  BearGPS$datetime <- strptime(paste(date,time),"%Y-%m-%d %H:%M:%S", tz="GMT")
  BearGPS$LAT <- unlist(as.vector(BearDATA[colnames[LAT_index]]),use.names = FALSE)
  BearGPS$LONG <- unlist(as.vector(BearDATA[colnames[LON_index]]),use.names = FALSE)
  if (length(ALT_index)>0) {
    BearGPS$ALT <- unlist(as.vector(BearDATA[colnames[ALT_index]]),use.names = FALSE)
  } else {
    BearGPS$ALT <- NA
  }
  BearGPS$LAT <- as.numeric(BearGPS$LAT)
  BearGPS$LONG <- as.numeric(BearGPS$LONG)
  BearGPS$ALT <- as.numeric(BearGPS$ALT)
  BearGPS <- BearGPS[order(BearGPS$datetime),]
  # remove ouliers:
  ALT <- BearGPS$ALT
  if(!all(is.na(ALT))) {
    z_scores <- abs(ALT-mean(ALT,na.rm=TRUE))/sd(ALT,na.rm=TRUE)
    BearGPS <- BearGPS[!(z_scores>3), ]
  } else {
    message("Warning: altitude missing")
  }
  BearGPS <- BearGPS[!is.na(BearGPS$LAT)&BearGPS$LAT>0,]
  LAT <- BearGPS$LAT
  z_scores <- abs(LAT-mean(LAT,na.rm=TRUE))/sd(LAT,na.rm=TRUE)
  BearGPS <- BearGPS[!(z_scores>5), ]
  LONG <- BearGPS$LONG
  z_scores <- abs(LONG-mean(LONG,na.rm=TRUE))/sd(LONG,na.rm=TRUE)
  BearGPS <- BearGPS[!(z_scores>5), ]
  filterfirst <- c(rep(TRUE,nrow(BearGPS)-1),FALSE)
  filterlast <- c(FALSE,rep(TRUE,nrow(BearGPS)-1))
  BearGPS$dist <- NA
  BearGPS$dtime <- NA
  BearGPS$dist[filterlast] <- distGeo(BearGPS[filterfirst,c("LONG","LAT")],BearGPS[filterlast,c("LONG","LAT")])
  BearGPS$dtime[filterlast] <- difftime(BearGPS$datetime[filterlast],BearGPS$datetime[filterfirst],units = "hours")
  BearGPS$climb <- NA
  if(!all(is.na(ALT))) {
    BearGPS$climb[filterlast] <- BearGPS$ALT[filterlast]-BearGPS$ALT[filterfirst]
  }
  
  return(BearGPS) 
}

BearGPS <- process_beardata(GPS_files[1])
BearGPS$track_ID <- 1
BearGPS$bear_name <- Bear_names[1]

for(i in 2:length(GPS_files)) {
  nextbear <- process_beardata(GPS_files[i])
  nextbear$track_ID <- i
  nextbear$bear_name <- Bear_names[i]
  
  BearGPS <- rbind(BearGPS,nextbear)
}

# make data spatial:
coordinates(BearGPS) <- ~LONG+LAT
proj4string(BearGPS) <- CRS("+init=epsg:4326")

# Move to new dataset:
ALL_Bears <- BearGPS
# remove small time intervals:
ALL_Bears <- ALL_Bears[ALL_Bears$dtime>0.1|is.na(ALL_Bears$dtime),]

Bear_subset <- ALL_Bears[ALL_Bears$bear_name=="Spiridon",]
#Bear_subset <- ALL_Bears[sample(1:nrow(ALL_Bears),5000),]

# map data:
palette(colorRampPalette(c("darkgreen","lightgreen","wheat","darkgrey"))(12))
ALT_pal <- colorFactor(palette(),domain=Bear_subset$ALT)

leaflet() %>% 
  addTiles(group="OSM (default)") %>% 
  addProviderTiles(providers$OpenTopoMap,group="Open topo map") %>% 
  addProviderTiles(providers$Esri.WorldTopoMap,group="ESRI topo map") %>%  
  addCircles(data=Bear_subset,color=~ALT_pal(Bear_subset$ALT)) %>%
  addMeasure(primaryLengthUnit = "metres")%>% 
  addLayersControl(baseGroups = c("OSM (default)","Open topo map","ESRI topo map"))

lags <- round(seq(1.4,5,by=0.4)^2)


###########
# Load DEM

# The DEM is provided by Copernicus and is the COP-DEM_GLO-30-DGED
# https://doi.org/10.5270/ESA-c5d3d65
Altitude_map <- raster("DEM_Stara_Planina.tif")
Altitude_map <- Altitude_map*10

Bear_subset$ALT_extract <- extract(Altitude_map,Bear_subset)
leaflet() %>% 
  addTiles(group="OSM (default)") %>% 
  addProviderTiles(providers$OpenTopoMap,group="Open topo map") %>% 
  addProviderTiles(providers$Esri.WorldTopoMap,group="ESRI topo map") %>%  
  addRasterImage(x=Altitude_map, 
                 col=palette(),
                 opacity = .7) %>% 
  addCircles(data=Bear_subset,color=~ALT_pal(Bear_subset$ALT)) %>%
  addMeasure(primaryLengthUnit = "metres")%>% 
  addLayersControl(baseGroups = c("OSM (default)","Open topo map","ESRI topo map"))

#########
# Combine Data
BearMove <- NULL

for(i in 1:length(GPS_files)) {
  BearGPS <- ALL_Bears[ALL_Bears$track_ID==i,]
  
  N <- nrow(BearGPS)
  
  interval <- 25    # for balanced sampling
  lag <- 1          # defining delta time
  filterfirst <- 1:N%%(interval) == 1
  filterlast <- 1:N%%(interval) == (lag+1)%%interval
  # if lag is larger than interval:
  Nremove <- lag%/%interval
  if(Nremove>0) {
    filterlast[filterlast][1:Nremove] <- FALSE
    filterfirst[filterfirst][(length(filterfirst[filterfirst])-Nremove+1):length(filterfirst[filterfirst])] <- FALSE
  }
  # if first and last are not equal:
  if(sum(filterfirst)!=sum(filterlast)){filterfirst[max(which(filterfirst))] <- FALSE}
  
  distBear <- distGeo(BearGPS[filterfirst,],BearGPS[filterlast,])
  endtime <- BearGPS$datetime[filterlast]
  starttime <- BearGPS$datetime[filterfirst]
  dtimeBear <- difftime(endtime,starttime,units = "hours")
  climbBear <- BearGPS$ALT[filterlast]-BearGPS$ALT[filterfirst]
  monthBear <- month(BearGPS$datetime)[filterfirst]
  
  #lag2:
  #lag <- 2
  for(lag in lags) {
    FF <- 1:N%%(interval) == 1
    FL <- 1:N%%(interval) == (lag+1)%%interval
    # if lag is larger than interval:
    Nremove <- lag%/%interval
    if(Nremove>0) {
      FL[FL][1:Nremove] <- FALSE
      FF[FF][(length(FF[FF])-Nremove+1):length(FF[FF])] <- FALSE
    }
    # if first and last are not equal:
    if(sum(FF)!=sum(FL)){FF[max(which(FF))] <- FALSE}
    dist <- distGeo(BearGPS[FF,],BearGPS[FL,])
    
    end <- BearGPS$datetime[FL]
    start <- BearGPS$datetime[FF]
    dtime <- difftime(end,start,units = "hours")
    climb <- BearGPS$ALT[FL]-BearGPS$ALT[FF]
    month <- month(BearGPS$datetime)[FF]
    distBear <- c(distBear,dist)
    dtimeBear <- c(dtimeBear,dtime)
    climbBear <- c(climbBear,climb)
    monthBear <- c(monthBear,month)
    starttime <- c(starttime,start)
    endtime <- c(endtime,end)
  }
  dtimeBear <- as.numeric(dtimeBear)
  temp <- data.frame(start_time=starttime,
                     end_time=endtime,
                     dtime=dtimeBear,
                     dist=distBear,
                     dist_km=distBear/1000,
                     climb=climbBear,
                     month=monthBear,
                     track_ID=i,
                     bear_name=Bear_names[i])
  # remove stationary observations - could result from detached collar or hibernation
  temp <- temp[temp$dist!=0,]
  
  BearMove <- rbind(BearMove,temp)
  remove(temp)
}
# remove NA's in dtime:
BearMove <- BearMove[!is.na(BearMove$dtime),]
BearMove <- BearMove[BearMove$dtime>0,]
BearMove$bear_name <- factor(BearMove$bear_name)
# male or female
Names <- levels(BearMove$bear_name)
Gender <- c("M","F","F","M","M","M","F","F","M","F","F","M","M","M","M","F","F","M","M","F")
MaleNames <- Names[Gender=="M"]
BearMove$Gender <- factor(ifelse(BearMove$bear_name %in% MaleNames,"Male","Female"))
BearMove <- BearMove[BearMove$month>3&BearMove$month<12,]
BearMove$log_dtime <- log(BearMove$dtime)

# colour per individual:
red <- data.frame(red=c(0,0.5,1))
green <- data.frame(green=c(0,0.5,1))
blue <- data.frame(blue=c(0,0.5,1))
coloursRG <- merge(red,green)
coloursRGB <- merge(coloursRG,blue)
palette(rgb(coloursRGB$red,coloursRGB$green,coloursRGB$blue,alpha=1))

plot(log(dist_km)~log(dtime/24),BearMove,
     #xlim=c(0,40),
     #ylim=c(0,15),
     pch=20, cex=0.5, col="grey70",
     main="Bear displacement",
     xlab="Days since last observation (log)", 
     ylab="Distance from origin (log(km))")
lines(aggregate(log(dist_km)~round(log(dtime/24),1),BearMove,mean),lwd=2)
lines(aggregate(log(dist_km)~round(log(dtime/24),1),BearMove,function(x){quantile(x,0.9772)}),lty=2, lwd=2)
#lines(aggregate(dist_km~round(dtime),BearMove,function(x){quantile(x,0.05)}),lty=2, lwd=2)
legend("topleft",legend=c("raw data","median per day","upper 95% bound per day"),
       lty=c(NA,1,2),lwd=c(NA,2,2),pch=c(20,NA,NA),col=c("grey70","black","black"))
Model_meanmove <- lm(log(dist_km)~log(dtime/24),BearMove)
abline(Model_meanmove,col="red")
sumstat <- aggregate(log(dist_km)~round(log(dtime/24),1)+Gender,BearMove,function(x){quantile(x,0.9772)})
sumstat$N <-  aggregate(log(dist_km)~round(log(dtime/24),1)+Gender,BearMove,length)[,3]
names(sumstat) <- c("log_dtime","Gender","log_dist","N")
sumstat <- sumstat[sumstat$log_dtime <4,]
#sumstat <- sumstat[sumstat$dtime >0,]
Model_uppermove <- lm(log_dist~log_dtime+Gender,sumstat,weights = sumstat$N)
summary(Model_uppermove)
abline(Model_uppermove,col="blue")

# make clusters to correct for bimodal distribution in movement data:
library(mclust)
palette(c("black","red"))
gmm_model <- Mclust(log(BearMove[c("dist_km")]), G = 2) 
bins <- predict(gmm_model)$classification
BearMove$cluster <- bins
# split in two:
subset <- sample(1:nrow(BearMove),5000)
plot(log(dist_km)~log(dtime/24),BearMove[subset,],
     #xlim=c(0,40),
     #ylim=c(0,15),
     pch=20, cex=0.5, col=BearMove$cluster[subset],
     main="Bear displacement",
     xlab="Hours since last observation (log)", 
     ylab="Distance from origin (log(km))")
Model_move1 <- mean(log(BearMove$dist_km[BearMove$cluster==1]))
Model_move2 <- lm(log(dist_km)~log(dtime/24),BearMove[BearMove$cluster==2,])
#Model_move2 <- mean(log(BearMove$dist_km)[BearMove$cluster==2])
abline(h=Model_move1,col="black")
abline(Model_move2,col="red")


bins_prior <- 0
bins_posterior <- 1
while(sum(bins_prior!=bins_posterior)>0) {
  bins_prior <- bins
  p1 <- dnorm(log(BearMove$dist_km)-Model_move1,sd=sd(log(BearMove$dist_km)-Model_move1))
  p2 <- dnorm(log(BearMove$dist_km)-predict(Model_move2,BearMove),sd=sd(log(BearMove$dist_km)-predict(Model_move2,BearMove),na.rm=TRUE))
  bins <- (p1 < p2) +1
  BearMove$cluster <- bins
  subset <- sample(1:nrow(BearMove),5000)
  plot(log(dist_km)~log(dtime/24),BearMove[subset,],
       #xlim=c(0,40),
       #ylim=c(0,15),
       pch=20, cex=0.5, col=BearMove$cluster[subset],
       main="Bear displacement",
       xlab="Hours since last observation (log)", 
       ylab="Distance from origin (log(km))")
  Model_move1 <- mean(log(BearMove$dist_km[BearMove$cluster==1]))
  Model_move2 <- lm(log(dist_km)~log(dtime/24),BearMove[BearMove$cluster==2,])
  #Model_move2 <- mean(log(BearMove$dist_km)[BearMove$cluster==2])
  abline(h=Model_move1,col="black")
  abline(Model_move2,col="red")
  bins_posterior <- bins
}
summary(as.factor(BearMove$cluster))

bins <- (order(BearMove$dtime[BearMove$cluster==2])-1)%/%250
sumstat <- aggregate(log(dtime)~bins,BearMove[BearMove$cluster==2,],median)
sumstat$log_dist <- aggregate(log(dist_km)~bins,BearMove[BearMove$cluster==2,],function(x){quantile(x,0.9772^2)})[,2]
sumstat$N <-  aggregate(log(dist_km)~bins,BearMove[BearMove$cluster==2,],length)[,2]
sumstat$median_dist <-  aggregate(log(dist_km)~bins,BearMove[BearMove$cluster==2,],median)[,2]
names(sumstat) <- c("bin","log_dtime","upper_dist","N","median_dist")
sumstat <- sumstat[sumstat$N>200,]    # remove long tail with insufficient data
sumstat$dtime <- exp(sumstat$log_dtime)
sumstat$time_dist_upp <- sumstat$dtime*exp(sumstat$upper_dist)
sumstat$time_dist_med <- sumstat$dtime*exp(sumstat$median_dist)
# fit model
Model_uppermove <- lm(upper_dist~log_dtime,sumstat,weights = sumstat$N)
Model_move2 <- lm(median_dist~log_dtime,sumstat,weights = sumstat$N)
summary(Model_uppermove)
summary(Model_move2)

# plot median and upper 95% bound
plot(dist_km~dtime,BearMove[BearMove$cluster==2,],xlim=c(0,10),ylim=c(0,7),pch=20, cex=0.5, col="grey70",main="Bear displacement",xlab="Hours since last observation", ylab="Distance from origin (km)")
lines(exp(sumstat$log_dtime),exp(sumstat$median_dist),lwd=2)
lines(exp(sumstat$log_dtime),exp(sumstat$upper_dist),lwd=2,lty=2)
#lines(aggregate(dist_km~round(dtime),BearMove,function(x){quantile(x,0.05)}),lty=2, lwd=2)
legend("topleft",legend=c("raw data","median per day","upper 95% bound per day","predicted"),cex=0.8,
       lty=c(NA,1,2,1),lwd=c(NA,2,2,1),pch=c(20,NA,NA,NA),col=c("grey70","black","black","blue"))
hours <- c(seq(0.01,0.99,by=0.01),1:150)
distance <- exp(coef(Model_move2)[1]+coef(Model_move2)[2]*log(hours))
lines(distance~hours,col="blue",lwd=2)
distanceF <- exp(coef(Model_uppermove)[1]+coef(Model_uppermove)[2]*log(hours))
lines(distanceF~hours,col="blue",lwd=2,lty=2)

# plot median and upper 95% bound
plot(dist_km~dtime,BearMove[BearMove$cluster==2,],xlim=c(0,144),ylim=c(0,15),pch=20, cex=0.5, col="grey70",main="Bear displacement",xlab="Hours since last observation", ylab="Distance from origin (km)")
legend("topleft",legend=c("raw data","median","upper 95% bound","predicted rational","predicted exponential"),cex=0.8,
       lty=c(NA,1,2,1,2),lwd=c(NA,2,2,2,2),pch=c(20,NA,NA,NA,NA),col=c("grey70","black","black","blue","blue"))
lines(exp(sumstat$log_dtime),exp(sumstat$median_dist),lwd=2)
lines(exp(sumstat$log_dtime),exp(sumstat$upper_dist),lwd=2,lty=2)
# try different model:
Model_uppermove2 <- lm(exp(upper_dist)~time_dist_upp+dtime+0,sumstat,weights = sumstat$N/exp(sumstat$upper_dist))
Model_meanmove2 <- lm(exp(median_dist)~time_dist_med+dtime+0,sumstat,weights = sumstat$N/exp(sumstat$upper_dist))
hours <- exp(seq(-8,5,by=0.1))
b <- -1/coef(Model_uppermove2)[1]
a <- coef(Model_uppermove2)[2]*b
distance <- a*hours/(b+hours)
lines(distance~hours,col="blue",lwd=2)
b <- -1/coef(Model_meanmove2)[1]
a <- coef(Model_meanmove2)[2]*b
distance <- a*hours/(b+hours)
lines(distance~hours,col="blue",lwd=2)
distance <- exp(coef(Model_move2)[1]+coef(Model_move2)[2]*log(hours))
lines(distance~hours,col="blue",lty=2,lwd=2)
distanceM <- exp(coef(Model_uppermove)[1]+coef(Model_uppermove)[2]*log(hours))
lines(distanceM~hours,col="blue",lty=2,lwd=2)

## Leave three bears out Cross-Validation: (= 7-fold CV)
TID <- unique(BearMove$track_ID)
folds <- (1:length(TID) - 1) %/% 3 + 1
set.seed(5000) # set random seed for fixed results
TID <- sample(TID)
CV_results <- data.frame()
for (FOLD in 1:7) {
  ID_test <- TID[folds==FOLD]
  ID_train <- TID[folds!=FOLD]
  
  BearMove.subset <- BearMove[BearMove$track_ID%in%ID_train&BearMove$cluster==2,]
  BearMove.subset <- BearMove.subset[order(BearMove.subset$dtime),]
  BearMove.valset <- BearMove[BearMove$track_ID%in%ID_test&BearMove$cluster==2,]
  BearMove.valset <- BearMove.valset[order(BearMove.valset$dtime),]
  
  # generate training set
  bins <- (order(BearMove.subset$dtime)-1)%/%250
  sumstat.train <- aggregate(log(dtime)~bins,BearMove.subset,median)
  sumstat.train$log_dist <- aggregate(log(dist_km)~bins,BearMove.subset,function(x){quantile(x,0.9772)})[,2]
  sumstat.train$N <-  aggregate(log(dist_km)~bins,BearMove.subset,length)[,2]
  sumstat.train$median_dist <-  aggregate(log(dist_km)~bins,BearMove.subset,median)[,2]
  names(sumstat.train) <- c("bin","log_dtime","upper_dist","N","median_dist")
  sumstat.train <- sumstat.train[sumstat.train$N>200,]    # remove long tail with insufficient data
  sumstat.train$dtime <- exp(sumstat.train$log_dtime)
  sumstat.train$time_dist_upp <- sumstat.train$dtime*exp(sumstat.train$upper_dist)
  sumstat.train$time_dist_med <- sumstat.train$dtime*exp(sumstat.train$median_dist)
  
  # fit model
  Model_uppermove <- lm(upper_dist~log_dtime,sumstat.train,weights = sumstat.train$N)
  Model_medianmove <- lm(median_dist~log_dtime,sumstat.train,weights = sumstat.train$N)
  Model_uppermove2 <- lm(exp(upper_dist)~time_dist_upp+dtime+0,sumstat.train,weights = sumstat.train$N/exp(sumstat.train$upper_dist))
  Model_meanmove2 <- lm(exp(median_dist)~time_dist_med+dtime+0,sumstat.train,weights = sumstat.train$N/exp(sumstat.train$upper_dist))
  
  # generate test set
  bins <- (order(BearMove.valset$dtime)-1)%/%100
  sumstat.test <- aggregate(log(dtime)~bins,BearMove.valset,median)
  sumstat.test$log_dist <- aggregate(log(dist_km)~bins,BearMove.valset,function(x){quantile(x,0.9772)})[,2]
  sumstat.test$N <-  aggregate(log(dist_km)~bins,BearMove.valset,length)[,2]
  sumstat.test$median_dist <-  aggregate(log(dist_km)~bins,BearMove.valset,median)[,2]
  names(sumstat.test) <- c("bin","log_dtime","upper_dist","N","median_dist")
  sumstat.test <- sumstat.test[sumstat.test$N>80,]    # remove long tail with insufficient data
  sumstat.test$dtime <- exp(sumstat.test$log_dtime)
  sumstat.test$time_dist_upp <- sumstat.test$dtime*exp(sumstat.test$upper_dist)
  sumstat.test$time_dist_med <- sumstat.test$dtime*exp(sumstat.test$median_dist)
  
  # add predicted values
  sumstat.test$upper_pred <- predict(Model_uppermove,sumstat.test)
  sumstat.test$median_pred <- predict(Model_medianmove,sumstat.test)
  b <- -1/coef(Model_uppermove2)[1]
  a <- coef(Model_uppermove2)[2]*b
  sumstat.test$upper_pred2 <- log(a*sumstat.test$dtime/(b+sumstat.test$dtime))
  b <- -1/coef(Model_meanmove2)[1]
  a <- coef(Model_meanmove2)[2]*b
  sumstat.test$median_pred2 <- log(a*sumstat.test$dtime/(b+sumstat.test$dtime))
  
  # calculate squared error and total variance: 
  sumstat.test$upper_se1 <- (sumstat.test$upper_dist-sumstat.test$upper_pred)^2
  sumstat.test$upper_se2 <- (sumstat.test$upper_dist-sumstat.test$upper_pred2)^2
  sumstat.test$upper_var <- (sumstat.test$upper_dist-mean(sumstat.test$upper_dist))^2
  sumstat.test$median_se1 <- (sumstat.test$median_dist-sumstat.test$median_pred)^2
  sumstat.test$median_se2 <- (sumstat.test$median_dist-sumstat.test$median_pred2)^2
  sumstat.test$median_var <- (sumstat.test$median_dist-mean(sumstat.test$median_dist))^2
  sumstat.test$fold <- FOLD
  
  CV_results <- rbind(CV_results,sumstat.test)
  
}

# calculate performance metrics
Rsquared_upper <- 1-sum(CV_results$upper_se1)/sum(CV_results$upper_var)
Rsquared_median <- 1-sum(CV_results$median_se1)/sum(CV_results$median_var)
Rsquared_upper_adj <- 1-(1-Rsquared_upper)*(nrow(CV_results)-1)/(nrow(CV_results)-2)
Rsquared_median_adj <- 1-(1-Rsquared_median)*(nrow(CV_results)-1)/(nrow(CV_results)-2)
Rsquared_upper2 <- 1-sum(CV_results$upper_se2)/sum(CV_results$upper_var)
Rsquared_median2 <- 1-sum(CV_results$median_se2)/sum(CV_results$median_var)
Rsquared_upper_adj2 <- 1-(1-Rsquared_upper2)*(nrow(CV_results)-1)/(nrow(CV_results)-2)
Rsquared_median_adj2 <- 1-(1-Rsquared_median2)*(nrow(CV_results)-1)/(nrow(CV_results)-2)
Rsquared_upper_adj
Rsquared_median_adj
Rsquared_upper_adj2
Rsquared_median_adj2

# Model for climbing - warning! not all altitude measurements are accurate
plot(log(abs(climb))~log(dtime),BearMove,
     #xlim=c(0,40),
     #ylim=c(0,15),
     pch=20, cex=0.5, col="grey70",
     main="Bear climb/descend",
     xlab="Hours since last observation (log)", 
     ylab="Altitude change (log(metres))")

bins <- (order(BearMove$dtime)-1)%/%250
legend("topleft",legend=c("raw data","upper 95% bound per day"),
       lty=c(NA,2),lwd=c(NA,2),pch=c(20,NA),col=c("grey70","black"))

sumstat <- aggregate(log(dtime)~bins,BearMove,median)
sumstat$log_climb <- aggregate(log(abs(climb))~bins,BearMove,function(x){quantile(x,0.9772)})[,2]
sumstat$N <-  aggregate(log(abs(climb))~bins,BearMove,length)[,2]
sumstat$median_climb <-  aggregate(log(abs(climb))~bins,BearMove,function(x){quantile(x,sqrt(.5))})[,2]
names(sumstat) <- c("bin","log_dtime","log_climb","N","median_climb")
sumstat <- sumstat[sumstat$N>200,]    # remove long tail with insufficient data
lines(median_climb~log_dtime,sumstat,lty=1, lwd=2)
lines(log_climb~log_dtime,sumstat,lty=2, lwd=2)
# fit model
Model_upperclimb <- lm(log_climb~log_dtime,sumstat,weights = sumstat$N)
Model_medianclimb <- lm(median_climb~log_dtime,sumstat,weights = sumstat$N)
summary(Model_upperclimb)
summary(Model_medianclimb)
abline(Model_upperclimb,col="blue")
abline(Model_medianclimb,col="red")

## Leave three bears out Cross-Validation: (= 6-fold CV)
TID <- unique(BearMove$track_ID[!is.na(BearMove$climb)])
folds <- (1:length(TID) - 1) %/% 3 + 1
TID <- sample(TID)
CV_results <- data.frame()
for (FOLD in 1:6) {
  ID_test <- TID[folds==FOLD]
  ID_train <- TID[folds!=FOLD]
  
  BearMove.subset <- BearMove[BearMove$track_ID%in%ID_train,]
  BearMove.subset <- BearMove.subset[order(BearMove.subset$dtime),]
  BearMove.valset <- BearMove[BearMove$track_ID%in%ID_test,]
  BearMove.valset <- BearMove.valset[order(BearMove.valset$dtime),]
  
  # generate training set
  bins <- (order(BearMove.subset$dtime)-1)%/%250
  sumstat.train <- aggregate(log(dtime)~bins,BearMove.subset,median)
  sumstat.train$log_climb <- aggregate(log(abs(climb))~bins,BearMove.subset,function(x){quantile(x,0.9772^2)})[,2]
  sumstat.train$N <-  aggregate(log(abs(climb))~bins,BearMove.subset,length)[,2]
  sumstat.train$median_climb <-  aggregate(log(abs(climb))~bins,BearMove.subset,function(x){quantile(x,sqrt(.5))})[,2]
  names(sumstat.train) <- c("bin","log_dtime","log_climb","N","median_climb")
  sumstat.train <- sumstat.train[sumstat.train$N>200,]    # remove long tail with insufficient data
  
  # fit model
  Model_upperclimb <- lm(log_climb~log_dtime,sumstat.train,weights = sumstat.train$N)
  Model_medianclimb <- lm(median_climb~log_dtime,sumstat.train,weights = sumstat.train$N)
  
  # generate test set
  bins <- (order(BearMove.valset$dtime)-1)%/%100
  sumstat.test <- aggregate(log(dtime)~bins,BearMove.valset,median)
  sumstat.test$log_climb <- aggregate(log(abs(climb))~bins,BearMove.valset,function(x){quantile(x,0.9772^2)})[,2]
  sumstat.test$N <-  aggregate(log(abs(climb))~bins,BearMove.valset,length)[,2]
  sumstat.test$median_climb <-  aggregate(log(abs(climb))~bins,BearMove.valset,function(x){quantile(x,sqrt(.5))})[,2]
  names(sumstat.test) <- c("bin","log_dtime","log_climb","N","median_climb")
  sumstat.test <- sumstat.test[sumstat.test$N>80,]    # remove long tail with insufficient data
  # add predicted values
  sumstat.test$upper_pred <- predict(Model_upperclimb,sumstat.test)
  sumstat.test$median_pred <- predict(Model_medianclimb,sumstat.test)
  # calculate squared error and total variance: 
  sumstat.test$upper_se <- (sumstat.test$log_climb-sumstat.test$upper_pred)^2
  sumstat.test$upper_var <- (sumstat.test$log_climb-mean(sumstat.test$log_climb))^2
  sumstat.test$median_se <- (sumstat.test$median_climb-sumstat.test$median_pred)^2
  sumstat.test$median_var <- (sumstat.test$median_climb-mean(sumstat.test$median_climb))^2
  sumstat.test$fold <- FOLD
  
  CV_results <- rbind(CV_results,sumstat.test)
  
}

# plot data
plot(log(abs(climb))~log(dtime),BearMove[BearMove$track_ID%in%ID_train,],
     #xlim=c(0,40),
     #ylim=c(0,15),
     pch=20, cex=0.5, col="grey70",
     main="Bear climb/descend",
     xlim=c(-2,5),
     ylim=c(0,7),
     xlab="Hours since last observation (log)", 
     ylab="Altitude change (log(metres))")
lines(log_climb~log_dtime, sumstat.train, lty=2,lwd=2)
lines(median_climb~log_dtime, sumstat.train, lty=1,lwd=2)
#lines(aggregate(dist_km~round(dtime),BearMove,function(x){quantile(x,0.05)}),lty=2, lwd=2)
legend("topleft",legend=c("raw data","upper 95% bound per day"),
       lty=c(NA,2),lwd=c(NA,2),pch=c(20,NA),col=c("grey70","black"))
abline(Model_upperclimb,col="blue")
abline(Model_medianclimb,col="red")

#plot results:
plot(log(abs(climb))~log(dtime),BearMove[BearMove$track_ID%in%ID_test,],
     #xlim=c(0,40),
     #ylim=c(0,15),
     pch=20, cex=0.5, col="grey70",
     main="Bear climb/descend",
     xlim=c(-2,5),
     ylim=c(0,7),
     xlab="Hours since last observation (log)", 
     ylab="Altitude change (log(metres))")
lines(log_climb~log_dtime, sumstat.test, lty=2,lwd=2)
lines(median_climb~log_dtime, sumstat.test, lty=1,lwd=2)
abline(Model_upperclimb,col="blue")
abline(Model_medianclimb,col="red")
lines(upper_pred~log_dtime, sumstat.test, lty=2,lwd=2)
lines(median_pred~log_dtime, sumstat.test, lty=2,lwd=2)

# calculate performance metrics
Rsquared_upper <- 1-sum(CV_results$upper_se)/sum(CV_results$upper_var)
Rsquared_median <- 1-sum(CV_results$median_se)/sum(CV_results$median_var)
Rsquared_upper_adj <- 1-(1-Rsquared_upper)*(nrow(CV_results)-1)/(nrow(CV_results)-2)
Rsquared_median_adj <- 1-(1-Rsquared_median)*(nrow(CV_results)-1)/(nrow(CV_results)-2)
Rsquared_upper_adj
Rsquared_median_adj

# check how many observations fall outside 95% confidence interval overall:
summary(abs(BearMove$climb)>exp(predict(Model_upperclimb,BearMove)))
3619/(3619+67650)

# plot with absolute numbers
plot(climb~dtime,BearMove,xlim=c(0,72),ylim=c(-750,750),pch=20, cex=0.5, col="grey70",main="Bear climb",xlab="Hours since last observation", ylab="Altitude change (metres)")
lines(exp(sumstat$log_dtime),exp(sumstat$log_climb),lty=2,lwd=2)
lines(exp(sumstat$log_dtime),-exp(sumstat$log_climb),lty=2,lwd=2)
lines(exp(sumstat$log_dtime),exp(sumstat$median_climb),lwd=2)
lines(exp(sumstat$log_dtime),-exp(sumstat$median_climb),lwd=2)
#lines(aggregate(dist_km~round(dtime),BearMove,function(x){quantile(x,0.05)}),lty=2, lwd=2)
legend("topleft",legend=c("raw data","50% bounds","95% bounds","predicted"),cex=0.8,
       lty=c(NA,1,2,1),lwd=c(NA,2,2,2),pch=c(20,NA,NA,NA),col=c("grey70","black","black","blue"))
days <- c(seq(0.01,0.99,by=0.01),1:150)
climbM <- exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(days))
lines(climbM~days,col="blue",lwd=2,lty=2)
climbM <- exp(coef(Model_medianclimb)[1]+coef(Model_medianclimb)[2]*log(days))
lines(climbM~days,col="blue",lwd=2)
climbM <- -exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(days))
lines(climbM~days,col="blue",lwd=2,lty=2)
climbM <- -exp(coef(Model_medianclimb)[1]+coef(Model_medianclimb)[2]*log(days))
lines(climbM~days,col="blue",lwd=2)

hours <- seq(.001,24,by=0.01)
dist_M2 <- exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(hours))
plot(dist_M2~hours, type="l")
abline(v=1,lty=2)
abline(h= exp(coef(Model_upperclimb)[1]),lty=2)
abline(v=5,lty=2)
abline(h= exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(5)),lty=2)
abline(v=1/60,lty=2)
abline(h= exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(1/60)),lty=2)


###########################
# Putting it all together #
###########################

# use Model_medianclimb, Model_upperclimb, Model_meanmove2 and Model_uppermove2
ALL_Bears$ALT_DEM <- extract(Altitude_map,ALL_Bears)
Bear_inarea <- ALL_Bears[!is.na(ALL_Bears$ALT_DEM)
                         &ALL_Bears$ALT_DEM!=253
                         &!is.na(ALL_Bears$datetime),]
Bear_inarea <- st_as_sf(Bear_inarea)
Bear_names <- unique(Bear_inarea$bear_name)

validation <- NULL

Bear <- Bear_names[1]
ITER <- 1
DT <- c(1,3,7,15,24,48)

for (DTIME in DT) {
  for (ITER in 1:50) {
    for (Bear in Bear_names) {
      VD <- data.frame(ITER=ITER,
                       Bear=Bear,
                       Dtime=NA,
                       Inradius_RED=NA,
                       Inradius_BLUE=NA,
                       Inpoly_RED=NA,
                       Inpoly_BLUE=NA,
                       RED_reduction=NA,
                       BLUE_reduction=NA)
      
      Bear_ind <- Bear_inarea[Bear_inarea$bear_name==Bear,]
      index <- FALSE
      while(!any(index)) {
        Bear_subset <- Bear_ind[sample(1:nrow(Bear_ind),1),]
        difftime <- difftime(Bear_ind$datetime,
                             Bear_subset$datetime,unit = "hours")
        index <- difftime>0&abs(DTIME-difftime)==min(abs(DTIME-difftime))
      }
      next_obs <- Bear_ind[index,][1,]
      dtime <- as.numeric(difftime(next_obs$datetime,
                                   Bear_subset$datetime,unit = "hours"))
      VD$Dtime <- dtime
      
      #plot(ALT_DEM~ALT,Bear_subset)
      #abline(a=0,b=0.1,col="red")
      #dtime <- exp(sample(0:400/100,nrow(Bear_subset)))
      #Bear_subset <- ALL_Bears[ALL_Bears$N==106&ALL_Bears$track_ID==8,]
      #Bear_subset <- ALL_Bears[ALL_Bears$N==4575&ALL_Bears$track_ID==8,]
      #Bear_subset <- ALL_Bears[ALL_Bears$N==1107&ALL_Bears$track_ID==4,]
      # validation:
      upperrange_poly <- calcActivityArea(Bear_subset,dtime=dtime,
                                          BearAlt=Bear_subset$ALT_DEM,
                                          bounds="upper",
                                          Altitude_map = Altitude_map)
      prefrange_poly <- calcActivityArea(Bear_subset,dtime=dtime,
                                         BearAlt=Bear_subset$ALT_DEM,
                                         Altitude_map = Altitude_map)
      
      proj4string(upperrange_poly) <- CRS("+init=epsg:4326")
      proj4string(prefrange_poly) <- CRS("+init=epsg:4326")
      b <- 1/0.04773
      a <- 0.60903*b
      alpha <- 1.876
      beta <- 0.504
      distanceU <- a*dtime/(b+dtime)
      distanceU[dtime<10] <- exp(alpha+beta*log(dtime/24))[dtime<10] # use exponential dispersal model in first 24 hours
      distanceU <- distanceU*1000
      b <- 1/0.035137
      a <- 0.14364*b
      alpha <- 0.5909
      beta <- 0.4945
      distanceM <- a*dtime/(b+dtime)
      distanceM[dtime<10] <- exp(alpha+beta*log(dtime/24))[dtime<10] # use exponential dispersal model in first 24 hours
      distanceM <- distanceM*1000
      
      distance <- as.numeric(st_distance(Bear_subset,next_obs)[1,1])
      VD$Inradius_BLUE <- distance < distanceU
      VD$Inradius_RED <- distance < distanceM
      
      # difference old vs new:
      VD$BLUE_reduction <- area(upperrange_poly)/sum(distanceU*distanceU*pi)
      VD$RED_reduction <- area(prefrange_poly)/sum(distanceM*distanceM*pi)
      
      upperrange_poly <- st_as_sf(upperrange_poly)
      prefrange_poly <- st_as_sf(prefrange_poly)
      VD$Inpoly_BLUE <- length(st_intersects(next_obs,upperrange_poly)[[1]])
      VD$Inpoly_RED <- length(st_intersects(next_obs,prefrange_poly)[[1]])
      validation <- rbind(validation,VD)
      
    }
  }
}

validation <- validation[validation$Dtime>0.8,]

leaflet() %>% 
  #addTiles(group="OSM (default)") %>% 
  addProviderTiles(providers$OpenTopoMap,group="Open topo map") %>% 
  addProviderTiles(providers$Esri.WorldTopoMap,group="ESRI topo map") %>%  
  #addPolygons(data=distU_alt_poly) %>% 
  addPolygons(data=upperrange_poly,color = "blue") %>% 
  addPolygons(data=prefrange_poly,color="red") %>% 
  #addRasterImage(x=dist_in_altrange*(dist_in_altrange<distanceU),opacity = .5) %>% 
  addMarkers(data=Bear_subset,icon=BearIcon) %>%
  addCircles(data=Bear_subset,
             group = "Current activity radius",
             radius = distanceU,
             fill = NA, 
             weight = 2) %>% 
  addCircles(data=Bear_subset,
             group = "Current activity radius",
             radius = distanceM,
             fill = NA, color="red", 
             weight = 2) %>% 
  addCircles(data=next_obs,col="magenta",radius = 10)%>% 
  addMeasure(primaryLengthUnit = "metres")%>% 
  addLayersControl(baseGroups = c("Open topo map","ESRI topo map"))
