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
sumstat <- aggregate(log(dist_km)~round(log(dtime/24),1)+Gender,BearMove[BearMove$cluster==2,],function(x){quantile(x,0.9772)})
sumstat$N <-  aggregate(log(dist_km)~round(log(dtime/24),1)+Gender,BearMove[BearMove$cluster==2,],length)[,3]
names(sumstat) <- c("log_dtime","Gender","log_dist","N")
sumstat$median_dist <- aggregate(log(dist_km)~round(log(dtime/24),1)+Gender,BearMove[BearMove$cluster==2,],function(x){quantile(x,sqrt(.5))})[,3]
sumstat <- sumstat[sumstat$log_dtime <4,]
#sumstat <- sumstat[sumstat$dtime >0,]
Model_uppermove <- lm(log_dist~log_dtime+Gender,sumstat,weights = sumstat$N)
Model_move2 <- lm(median_dist~log_dtime+Gender,sumstat,weights = sumstat$N)
summary(Model_uppermove)
abline(Model_uppermove,col="blue")

# plot median and upper 95% bound
plot(dist_km~dtime,BearMove[BearMove$cluster==2,],xlim=c(0,144),ylim=c(0,15),pch=20, cex=0.5, col="grey70",main="Bear displacement",xlab="Hours since last observation", ylab="Distance from origin (km)")
Nobs <- aggregate(dist_km~round(log(dtime),1),BearMove[BearMove$cluster==2,],length)$dist_km
aggregated_median <- aggregate(dist_km~round(log(dtime/24),1),BearMove[BearMove$cluster==2,],function(x){quantile(x,sqrt(.5))})
names(aggregated_median) <- c("log_dtime","dist_km")
aggregated_median$dtime <- exp(aggregated_median$log_dtime)*24
lines(dist_km~dtime,aggregated_median,lwd=2)
aggregated_upper <- aggregate(dist_km~round(log(dtime/24),1),BearMove[BearMove$cluster==2,],function(x){quantile(x,0.9772)})
names(aggregated_upper) <- c("log_dtime","dist_km")
aggregated_upper$dtime <- exp(aggregated_upper$log_dtime)*24
lines(dist_km~dtime,aggregated_upper,lty=2, lwd=2)
#lines(aggregate(dist_km~round(dtime),BearMove,function(x){quantile(x,0.05)}),lty=2, lwd=2)
legend("topleft",legend=c("raw data","median per day","upper 95% bound per day","predicted"),cex=0.8,
       lty=c(NA,1,2,1),lwd=c(NA,2,2,1),pch=c(20,NA,NA,NA),col=c("grey70","black","black","blue"))
hours <- c(seq(0.01,0.99,by=0.01),1:150)
distance <- exp(coef(Model_move2)[1]+coef(Model_move2)[2]*log(hours/24))
lines(distance~hours,col="darkgreen",lty=2)
distanceF <- exp(coef(Model_uppermove)[1]+coef(Model_uppermove)[2]*log(hours/24))
lines(distanceF~hours,col="magenta")
distanceM <- exp(coef(Model_uppermove)[1]+coef(Model_uppermove)[2]*log(hours/24)+coef(Model_uppermove)[3])
lines(distanceM~hours,col="blue")

# check how many observations fall outside 95% confidence interval overall:
summary(BearMove$dist_km[BearMove$month>3&BearMove$month<12]>exp(predict(Model_uppermove,BearMove)))
2213/(2213+76224)
2258/(2258+81825)


# Model for climbing - warning! not all altitude measurements are accurate
palette(rgb(coloursRGB$red,coloursRGB$green,coloursRGB$blue,alpha=1))

plot(log(abs(climb))~log(dtime),BearMove,
     #xlim=c(0,40),
     #ylim=c(0,15),
     pch=20, cex=0.5, col=BearMove$bear_name[BearMove$month>3&BearMove$month<12],
     main="Bear climb/descend",
     xlab="Hours since last observation (log)", 
     ylab="Altitude change (log(metres))")
lines(aggregate(log(abs(climb))~(round(log(dtime),1)),BearMove,function(x){quantile(x,sqrt(.5))}),lty=1, lwd=2)
lines(aggregate(log(abs(climb))~(round(log(dtime),1)),BearMove,function(x){quantile(x,0.9772^2)}),lty=2, lwd=2)
#lines(aggregate(dist_km~round(dtime),BearMove,function(x){quantile(x,0.05)}),lty=2, lwd=2)
legend("topleft",legend=c("raw data","upper 95% bound per day"),
       lty=c(NA,2),lwd=c(NA,2),pch=c(20,NA),col=c("grey70","black"))
sumstat <- aggregate(log(abs(climb))~round(log(dtime),1)+Gender,BearMove,function(x){quantile(x,0.9772^2)})
sumstat$N <-  aggregate(log(abs(climb))~round(log(dtime),1)+Gender,BearMove,length)[,3]
sumstat$median_climb <-  aggregate(log(abs(climb))~round(log(dtime),1)+Gender,BearMove,function(x){quantile(x,sqrt(.5))})[,3]
names(sumstat) <- c("log_dtime","Gender","log_climb","N","median_climb")
sumstat <- sumstat[sumstat$log_dtime <5,]
#sumstat <- sumstat[sumstat$dtime >0,]
Model_upperclimb <- lm(log_climb~log_dtime+Gender,sumstat,weights = sumstat$N)
Model_medianclimb <- lm(median_climb~log_dtime+Gender,sumstat,weights = sumstat$N)
summary(Model_upperclimb)
summary(Model_medianclimb)
abline(Model_upperclimb,col="blue")
abline(Model_medianclimb,col="red")

# check how many observations fall outside 95% confidence interval overall:
summary(abs(BearMove$climb)>exp(predict(Model_upperclimb,BearMove)))
3619/(3619+67650)

# plot median and upper 95% bound
plot(climb~dtime,BearMove,xlim=c(0,144),ylim=c(-750,750),pch=20, cex=0.5, col="grey70",main="Bear climb",xlab="Hours since last observation", ylab="Altitude change (metres)")
Nobs <- aggregate(climb~round(log(dtime),1),BearMove,length)$climb
aggregated_upper <- aggregate(climb~round(log(dtime),1),BearMove,function(x){quantile(x,0.9772)})[Nobs>100,]
names(aggregated_upper) <- c("log_dtime","climb")
aggregated_upper$dtime <- exp(aggregated_upper$log_dtime)
lines(climb~dtime,aggregated_upper,lty=2,lwd=2)
aggregated_upper$descend <- aggregate(climb~round(log(dtime),1),BearMove,function(x){quantile(x,1-0.9772)})[Nobs>100,]$climb
aggregated_upper$upper50 <- aggregate(climb~round(log(dtime),1),BearMove,function(x){quantile(x,sqrt(sqrt(.5)))})[Nobs>100,]$climb
aggregated_upper$lower50 <- aggregate(climb~round(log(dtime),1),BearMove,function(x){quantile(x,1-sqrt(sqrt(.5)))})[Nobs>100,]$climb
lines(descend~dtime,aggregated_upper,lty=2,lwd=2)
lines(lower50~dtime,aggregated_upper,lwd=2)
lines(upper50~dtime,aggregated_upper,lwd=2)
#lines(aggregate(dist_km~round(dtime),BearMove,function(x){quantile(x,0.05)}),lty=2, lwd=2)
legend("topleft",legend=c("raw data","upper 95% bound per day","predicted (male)","predicted (female)"),cex=0.8,
       lty=c(NA,2,1,1),lwd=c(NA,2,1,1),pch=c(20,NA,NA,NA),col=c("grey70","black","blue","magenta"))
days <- c(seq(0.01,0.99,by=0.01),1:150)
climbF <- exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(days))
lines(climbF~days,col="magenta")
climbF <- exp(coef(Model_medianclimb)[1]+coef(Model_medianclimb)[2]*log(days))
lines(climbF~days,col="magenta",lty=2)
climbM <- exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(days)+coef(Model_upperclimb)[3])
lines(climbM~days,col="blue")
climbM <- exp(coef(Model_medianclimb)[1]+coef(Model_medianclimb)[2]*log(days)+coef(Model_medianclimb)[3])
lines(climbM~days,col="blue",lty=2)
climbF <- -exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(days))
lines(climbF~days,col="magenta")
climbF <- -exp(coef(Model_medianclimb)[1]+coef(Model_medianclimb)[2]*log(days))
lines(climbF~days,col="magenta",lty=2)
climbM <- -exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(days)+coef(Model_upperclimb)[3])
lines(climbM~days,col="blue")
climbM <- -exp(coef(Model_medianclimb)[1]+coef(Model_medianclimb)[2]*log(days)+coef(Model_medianclimb)[3])
lines(climbM~days,col="blue",lty=2)


hours <- seq(.001,24,by=0.01)
dist_M2 <- exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(hours))
plot(dist_M2~hours, type="l")
abline(v=1,lty=2)
abline(h= exp(coef(Model_upperclimb)[1]),lty=2)
abline(v=5,lty=2)
abline(h= exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(5)),lty=2)
abline(v=1/60,lty=2)
abline(h= exp(coef(Model_upperclimb)[1]+coef(Model_upperclimb)[2]*log(1/60)),lty=2)

# try different model:
BearMove$time_dist <- BearMove$dtime*BearMove$dist_km
# plot median and upper 95% bound
plot(dist_km~dtime,BearMove[BearMove$cluster==2,],xlim=c(0,144),ylim=c(0,15),pch=20, cex=0.5, col="grey70",main="Bear displacement",xlab="Hours since last observation", ylab="Distance from origin (km)")
Nobs <- aggregate(dist_km~round(log(dtime),1),BearMove[BearMove$cluster==2,],length)$dist_km
aggregated_median <- aggregate(dist_km~round(log(dtime/24),1),BearMove[BearMove$cluster==2,],function(x){quantile(x,sqrt(.5))})
names(aggregated_median) <- c("log_dtime","dist_km")
aggregated_median$dtime <- exp(aggregated_median$log_dtime)*24
lines(dist_km~dtime,aggregated_median,lwd=2)
aggregated_upper <- aggregate(dist_km~round(log(dtime/24),1),BearMove[BearMove$cluster==2,],function(x){quantile(x,0.9772)})
names(aggregated_upper) <- c("log_dtime","dist_km")
aggregated_upper$dtime <- exp(aggregated_upper$log_dtime)*24
lines(dist_km~dtime,aggregated_upper,lty=2, lwd=2)
#lines(aggregate(dist_km~round(dtime),BearMove,function(x){quantile(x,0.05)}),lty=2, lwd=2)
legend("topleft",legend=c("raw data","median per day","upper 95% bound per day","predicted rational","predicted exponential"),cex=0.8,
       lty=c(NA,1,2,1,2),lwd=c(NA,2,2,2,2),pch=c(20,NA,NA,NA,NA),col=c("grey70","black","black","blue","blue"))
sumstat <- aggregate(dist_km~round(log(dtime),1),BearMove[BearMove$cluster==2,],function(x){quantile(x,0.9772)})
sumstat$N <-  aggregate(log(dist_km)~round(log(dtime),1),BearMove[BearMove$cluster==2,],length)[,2]
names(sumstat) <- c("log_dtime","dist","N")
sumstat$dtime <- exp(sumstat$log_dtime)
sumstat <- sumstat[sumstat$N>99,]
sumstat$time_dist <- sumstat$dtime*sumstat$dist
Model_uppermove2 <- lm(dist~time_dist+dtime+0,sumstat,weights = sumstat$N/sumstat$dist)
sumstat <- aggregate(dist_km~round(log(dtime),1),BearMove[BearMove$cluster==2,],function(x){quantile(x,sqrt(0.5))})
sumstat$N <-  aggregate(log(dist_km)~round(log(dtime),1),BearMove[BearMove$cluster==2,],length)[,2]
names(sumstat) <- c("log_dtime","dist","N")
sumstat$dtime <- exp(sumstat$log_dtime)
sumstat <- sumstat[sumstat$N>99,]
sumstat$time_dist <- sumstat$dtime*sumstat$dist
Model_meanmove2 <- lm(dist~time_dist+dtime+0,sumstat,weights = sumstat$N/sumstat$dist)
hours <- exp(seq(-8,5,by=0.1))
b <- -1/coef(Model_uppermove2)[1]
a <- coef(Model_uppermove2)[2]*b
distance <- a*hours/(b+hours)
lines(distance~hours,col="blue",lwd=2)
b <- -1/coef(Model_meanmove2)[1]
a <- coef(Model_meanmove2)[2]*b
distance <- a*hours/(b+hours)
lines(distance~hours,col="blue",lwd=2)
distance <- exp(coef(Model_move2)[1]+coef(Model_move2)[2]*log(hours/24))
lines(distance~hours,col="blue",lty=2,lwd=2)
distanceM <- exp(coef(Model_uppermove)[1]+coef(Model_uppermove)[2]*log(hours/24)+coef(Model_uppermove)[3])
lines(distanceM~hours,col="blue",lty=2,lwd=2)
for(sim in 1:50) {
  coefs <- summary(Model_uppermove2)$coef
  b <- -1/rnorm(1,coefs[1,1],coefs[1,2])
  a <- rnorm(1,coefs[2,1],coefs[2,2])*b
  distance <- a*hours/(b+hours)
  lines(distance~hours,col="magenta")
}

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
      b <- 1/0.04474
      a <- 0.71022*b
      alpha <- 2.0892
      beta <- 0.4872
      distanceU <- a*dtime/(b+dtime)
      distanceU[dtime<24] <- exp(alpha+beta*log(dtime/24))[dtime<24] # use exponential dispersal model in first 24 hours
      distanceU <- distanceU*1000
      b <- 1/0.05528
      a <- 0.30425*b
      alpha <- 1.08046
      beta <- 0.48566
      distanceM <- a*dtime/(b+dtime)
      distanceM[dtime<24] <- exp(alpha+beta*log(dtime/24))[dtime<24] # use exponential dispersal model in first 24 hours
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
