#Test GPS on a meerkat

library(dismo)
library(viridis)
library(lubridate)

source('~/Dropbox/meerkats/code/meerkat_functions.R')
source('~/Dropbox/Hyena_pilot/Code/make_google_earth_movies.R')

#directory of collar data
dir <- '/Users/astrandb/Desktop/Kalahari/HM_2017_2/COLLAR2/GPS2'

#output directory
outdir <- '/Users/astrandb/Desktop/Kalahari/Animations/anim_2017-08-25'

#filenames of collar data
filenames <- c('/HM_RT_A13_GPS_20170821-20170825/RT_A13_GPS_20170821-20170825.csv','/HM_LT_A05_GPS_20170821-20170825/HM_LT_A05_GPS_20170821-20170825.csv','/HM_HTB_A06_GPS_20170821-20170825/HM_HTB_A06_20170821-20170825.csv','/HM_HRT_A07_GPS_20170821-20170825/HM_HRT_A07_GPS_20170821-20170825.csv','/HM_HMB_A09_GPS_20170821-20170825/HM_HMB_A09_GPS_20170821-20170825.csv')

#dominant female file
df.file <- '/Users/astrandb/Desktop/Kalahari/HM_2017_2/HM_20170825/HM_VLF206_GPS_20170825/HM_VLF206_GPS_20170825.csv'

#other files
other.files <- c('/Users/astrandb/Desktop/Kalahari/HM_2017_2/HM_20170825/HM_HLT_GPS_20170825/HM_HLT_GPS_20170825.csv')

#date to plot
date.to.plot <- '2017-08-25'

#names of the meerkats, in order
meerkat.names <- c('RT','LT','HTB','HRT','PET','AJB','HLT')

#hours to plot
hours.to.plot <- c(8,9)

#list of all files
filenames.full <- c()
for(i in 1:length(filenames)){
	filenames.full <- c(filenames.full,paste(dir,filenames[i],sep=''))
}

filenames.full <- c(filenames.full,df.file,other.files)

gps.all.inds <- data.frame()
for(i in 1:length(filenames.full)){

	filename <- filenames.full[i]
	
	#load gps and clean up
	gps.all <- read.table(filename,header=TRUE,sep='\t')
	gps.all$datetime <- sapply(gps.all$timestamp,function(x){return(gsub('\\.00','',as.character(x)))})
	gps.all$datetime <- as.POSIXct(strptime(gps.all$datetime,'%d/%m/%Y %H:%M:%S'))
	gps.clean <- gps.all[which(!is.na(gps.all$datetime)),]
	
	#convert to utm
	utm <- latlon.to.utm(cbind(gps.clean$location.lon,gps.clean$location.lat),southern_hemisphere=T,utm.zone=32)
	gps.clean$easting <- utm[,1]
	gps.clean$northing <- utm[,2]

	gps.clean <- gps.clean[which(gps.clean$easting > 0),]
	gps.clean <- gps.clean[which(gps.clean$northing > 6900000),]
	gps.clean <- gps.clean[which(year(gps.clean$datetime)>2016),]
	
	gps.clean$id <- meerkat.names[i]
	
	gps.all.inds <- rbind(gps.all.inds,gps.clean)
}

#get relevant data for plotting
gps <- gps.all.inds[which(date(gps.all.inds$datetime)==date.to.plot & hour(gps.all.inds$datetime) %in% hours.to.plot),]

#get data into matrices for plotting
t0 <- min(gps$datetime)
tf <- max(gps$datetime)
times <- seq(t0,tf,1)
lats <- lons <- matrix(NA,nrow=length(meerkat.names),ncol=length(times))
for(i in 1:length(meerkat.names)){
	curr <- gps[which(gps$id==meerkat.names[i]),]
	lats[i,match(curr$datetime,times)] <- curr$location.lat
	lons[i,match(curr$datetime,times)] <- curr$location.lon
}

#plot trajectories on google earth
trajectories.movie(lats,lons,1,length(times),step=5,tail.time=120,base.dir='/Users/astrandb/Desktop',colors=rainbow(length(meerkat.names)),on.map=FALSE,ind.names=meerkat.names,plot.legend=T,times=times,show.all.past=FALSE,show.scale.bar=FALSE,scale.bar.len=50,utm.zone='34',southern_hemisphere=T,scale.bar.text='50 m',scale.bar.text.offset=NULL,playback.time=NULL,playback.lat=NULL,playback.lon=NULL,inds=NULL)


#set up plot
cols <- rainbow(length(meerkat.names))
symbols <- 18:(18+length(meerkat.names))
times <- unique(gps$datetime)
for(i in 1:length(times)){
	t <- times[i]
	png(filename=paste(outdir,'/',i,'.png',sep=''),height=800,width=800)
	plot(NULL,xlim=c(quantile(gps$easting,0.005,na.rm=T),quantile(gps$easting,0.995,na.rm=T)),ylim=c(quantile(gps$northing,0.005,na.rm=T),quantile(gps$northing,0.995,na.rm=T)),asp=1,xlab='Easting (m)',ylab='Northing (m)')
	legend('topright',legend=meerkat.names,fill=cols,col=cols,pch=symbols[j])
	gps.curr <- gps[which(as.character(gps$datetime)==t),]
	for(j in 1:length(meerkat.names)){
		id <- meerkat.names[j]
		gps.curr2 <- gps.curr[which(gps.curr$id==id),]
		points(gps.curr2$easting,gps.curr2$northing,col=cols[j],pch=symbols[j],cex=2)
	}
	dev.off()
	
}

plot(gps.clean$easting,gps.clean$northing,col=rainbow(nrow(gps.clean)),pch=19,cex=0.1,asp=1,xlab='Easting (m)',ylab='Northing (m)')

#Hour
t0 <- 1
dur <- 60*3
tf <- t0 + dur*60 
plot(gps.clean$easting[t0:tf],gps.clean$northing[t0:tf],type='l',xlab='Easting (m)',ylab='Northing (m)',asp=1)
points(gps.clean$easting[t0:tf],gps.clean$northing[t0:tf],pch=19,cex=0.3,col=rainbow(tf-t0+1))