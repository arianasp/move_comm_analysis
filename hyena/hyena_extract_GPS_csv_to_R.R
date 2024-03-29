#Extract hyena GPS data and save as R file

#cc16_352a: WRTH
#cc16_352b: BORA
#cc16_354a: BYTE
#cc16_360a: MGTA
#cc16_366a: FAY 

#whether to overwrite saved files on the server
overwrite <- F

#specify time zone offset from GMT
tz_offset <- 3 #3 hours

#Libraries
library(lubridate)

#Get useful functions
codedir <- '~/Dropbox/code_ari/move_comm_analysis/hyena'
setwd(codedir)
source('hyena_functions.R')

#Directory where csv files are stored
indir <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/rawdata/gps_technosmart'

#Directory where to store output
outdir <- '/Volumes/EAS_shared/hyena/archive/hyena_pilot_2017/processed/gps'

#Set working directory
setwd(indir)

#First create a data frame (hyena.ids) that maps hyenas to collar ids
hyena.ids <- data.frame(name = c('WRTH','BORA','BYTE','MGTA','FAY'), collar = c('352a','352b','354a','360a','366a'), id = seq(1,5,1),color=c('red','blue','green','magenta','orange'))
#Get all files in folder
files <- list.files()
files <- files[grep('.csv',files)]
hyena.ids$name <- as.character(hyena.ids$name)
hyena.ids$collar <- as.character(hyena.ids$collar)
hyena.ids$color <- as.character(hyena.ids$name)

#Read everything into a single data table (hyena.gps)
hyena.gps.list <- list()
for(i in 1:nrow(hyena.ids)){
	print(i)
	file <- paste('cc16_',hyena.ids$collar[i],'GPS.csv',sep='')
	collar <- hyena.ids$collar[i]
	dat.curr <- read.table(files[i],sep='\t',header=TRUE,stringsAsFactors=F)
	dat.curr$id <- i
	dat.curr$name <- hyena.ids$name[i]
	hyena.gps.list[[i]] <- dat.curr
}
hyena.gps <- do.call(rbind,hyena.gps.list)
rm(hyena.gps.list)

#add a datetime column
datetimes <- gsub('[.]00','',hyena.gps$timestamp)
hyena.gps$datetime <- strptime(datetimes,format='%d/%m/%Y %H:%M:%S',tz='UTC')

#Filter out data with NA timestamps
hyena.gps <- hyena.gps[which(!is.na(hyena.gps$datetime)),]

#Convert to UTM
utms <- latlon.to.utm(cbind(hyena.gps$location.long,hyena.gps$location.lat),utm.zone=36,southern_hemisphere=TRUE)
hyena.gps$easting <- utms[,1]
hyena.gps$northing <- utms[,2]

#Get list of timestamps (from first to last)
first.time <- floor_date(min(hyena.gps$datetime + tz_offset*60*60,na.rm=T), 'day') - tz_offset*60*60
last.time <- ceiling_date(max(hyena.gps$datetime + tz_offset*60*60 ,na.rm=T), 'day')  - tz_offset*60*60 - 1
timestamps <- seq(first.time,last.time,by=1,tz='UTC')

#Get matrices of lats and lons
n <- nrow(hyena.ids)
lats <- lons <- xs <- ys <- matrix(NA,nrow=n,ncol=length(timestamps))
for(i in 1:n){
	curr.dat <- hyena.gps[which(hyena.gps$id==i),]
	if(nrow(curr.dat)>0){
		idxs <- match(as.character(curr.dat$datetime),as.character(timestamps))
		non.nas <- which(!is.na(idxs))
		lats[i,idxs[non.nas]] <- curr.dat$location.lat[non.nas]
		lons[i,idxs[non.nas]] <- curr.dat$location.long[non.nas]
		xs[i,idxs[non.nas]] <- curr.dat$easting[non.nas]
		ys[i,idxs[non.nas]] <- curr.dat$northing[non.nas]
	}
}

#Get list of dates and indexes to starts of (complete) days (so skip the first Dec 31 hours)
days <- as.Date(timestamps + tz_offset*60*60)
dates <- sort(unique(days))
day.start.idxs <- c()
for(i in 1:length(dates)){
	day.start.idxs <- c(day.start.idxs, min(which(days==dates[i])))
}

if(overwrite){
  setwd(outdir)
  save(file='hyena_gps_level0.RData',list=c('xs','ys','lats','lons','timestamps','hyena.ids','hyena.gps','dates','day.start.idxs'))
  save(file='hyena_xy_level0.Rdata',list=c('xs','ys'))
  save(file='hyena_timestamps.Rdata',list=c('timestamps'))
  save(file='hyena_day_start_idxs.RData', list = c('day.start.idxs'))
  save(file='hyena_latlon_level0.RData',list=c('lats','lons'))
  save(file='hyena_ids.RData',list=c('hyena.ids'))
}
