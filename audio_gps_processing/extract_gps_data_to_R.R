#Read in GPS data from raw CSV files to an R data frame

#DIRECTORIES
collar.dirs <- c('~/Dropbox/meerkats/data/Kalahari2017/HM_2017/COLLAR/GPS','~/Dropbox/meerkats/data/Kalahari2017/HM_2017_2/COLLAR2/GPS2','~/Dropbox/meerkats/data/Kalahari2017/HM_2017_3/COLLAR3/GPS3')
focal.dirs <- c('~/Dropbox/meerkats/data/Kalahari2017/HM_2017','~/Dropbox/meerkats/data/Kalahari2017/HM_2017_2','~/Dropbox/meerkats/data/Kalahari2017/HM_2017_3')
garmin.files <- c('~/Dropbox/meerkats/data/Kalahari2017/HM_2017_2/HM_20170823/HM_HLT_GPS_20170823.txt','~/Dropbox/meerkats/data/Kalahari2017/HM_2017_2/HM_20170824/HM_HLT_GPS_20170824.txt')

#FUNCTIONS
source('~/Dropbox/meerkats/code/general/meerkat_functions.R')

get.start.date <- function(x){
	splitstr <- strsplit(x,'_')[[1]]
	dates <- splitstr[length(splitstr)]
	first.date <- strsplit(dates,'-')[[1]][1]
	return(first.date)
}
#FUNCTIONS
get.last.date <- function(x){
	splitstr <- strsplit(x,'_')[[1]]
	dates <- splitstr[length(splitstr)]
	last.date <- strsplit(dates,'-')[[1]][2]
	return(last.date)
}

#MAIN
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/ids.RData')

#Read in collar data
collar.data <- data.frame()
for(dir in collar.dirs){

	print(dir)
	#get subdirectories with GPS data
	subdirs <- list.files(dir)
	subdirs <- subdirs[grep('HM',subdirs)]
	
	#read in the csv file with the data
	for(i in 1:length(subdirs)){
		print(i)
		filename <- subdirs[i]
		path <- paste(dir,'/',filename,'/',filename,'.csv',sep='')
		dat <- read.table(path,header=T,sep='\t',stringsAsFactors=F)
		curr <- data.frame(file=rep(filename,nrow(dat)),timestamp=dat$timestamp,lat=dat$location.lat,lon=dat$location.lon)
		collar.data <- rbind(collar.data,curr)
	}
}

#convert to characters and numeric
collar.data$file <- as.character(collar.data$file)

#add a column for datetime
collar.data$datetime <- as.POSIXct(collar.data$timestamp,format='%d/%m/%Y %H:%M:%S')

#Get rid of intervals where collars were sitting in a shrub
collar.data$date <- as.Date(collar.data$datetime)
collar.data$date0 <- sapply(collar.data$file,get.start.date)
collar.data$datef <- sapply(collar.data$file,get.last.date)
collar.data$date0 <- as.Date(collar.data$date0,format='%Y%m%d')
collar.data$datef <- as.Date(collar.data$datef,format='%Y%m%d')
good.rows <- which(collar.data$date >= collar.data$date0 & collar.data$date <= collar.data$datef)
collar.data.good <- collar.data[good.rows,c('file','timestamp','lat','lon','datetime')]
collar.data.good$device <- 'collar'

#Read in focal gps data
focal.data <- data.frame()
for(dir in focal.dirs){
	print(dir)
	subdirs <- list.files(dir)
	subdirs <- subdirs[grep('HM',subdirs)]
	for(subdir in subdirs){
		subsubdirs <- list.files(paste(dir,'/',subdir,sep=''))
		to.remove <- unique(c(grep('.txt',subsubdirs),grep('SOUNDFOC',subsubdirs)))
		subsubdirs <- subsubdirs[-to.remove]
		for(subsubdir in subsubdirs){
			files <- list.files(paste(dir,'/',subdir,'/',subsubdir,sep=''))
			files <- files[grep('.csv',files)]
			for(file in files){
				filename <- file
				fullpath <- paste(dir,'/',subdir,'/',subsubdir,'/',file,sep='')
				dat <- read.table(file=fullpath,sep='\t',header=T,stringsAsFactors=F)
				curr <- data.frame(file=rep(filename,nrow(dat)),timestamp=dat$timestamp,lat=dat$location.lat,lon=dat$location.lon)
				focal.data <- rbind(focal.data,curr)
			}
		}
	}
}

focal.data$file <- as.character(focal.data$file)
focal.data$datetime <- as.POSIXct(focal.data$timestamp,format='%d/%m/%Y %H:%M:%S')
focal.data$device <- 'fluffy'

#get rid of focal data whose filename date does not match its datetime (i.e. stuff left on the tag from before)
datetime.str <- as.character(focal.data$datetime)
date.str <- sapply(datetime.str,FUN=function(x){return(strsplit(x,' ')[[1]][1])})
date.str.nodash <- gsub('-','',date.str)
date.str.both <- cbind(date.str.nodash,focal.data$file)
dates.match <- apply(date.str.both,1,FUN=function(x){return(grepl(x[1],x[2]))})
idxs.keep <- which(dates.match==1)

gps.data <- rbind(collar.data.good,focal.data[idxs.keep,])


#Get id and dye mark of focal
gps.data$id <- gps.data$dyemark <- gps.data$name <- NA
for(i in 1:nrow(meerkat.ids)){
	id <- as.character(meerkat.ids$id[i])
	dyemark <- as.character(meerkat.ids$dyemark[i])
	dyemark <- sub('\\+','',dyemark)
	name <- as.character(meerkat.ids$name[i])
	matches <- unique(c(grep(paste('_',id,'_',sep=''),gps.data$file,ignore.case=T),grep(paste('_',dyemark,'_',sep=''),gps.data$file,ignore.case=T),grep(paste('_',name,'_',sep=''),gps.data$file,ignore.case=T)))
	if(length(matches)>0){
		gps.data$id[matches] <- id
		gps.data$dyemark[matches] <- dyemark
		gps.data$name[matches] <- name
	}
}

#ADD DATA FROM GARMIN
garmin.data <- data.frame()
for(file in garmin.files){
	dat <- read.table(file,sep='\t',stringsAsFactors=F,header=F,fileEncoding="latin1")
	curr <- dat[,c(3,2)]
	colnames(curr) <- c('timestamp','pos')
	pos <- curr$pos
	curr$lat <- sapply(pos,FUN=function(x){return(-as.numeric(gsub('S','',strsplit(x,' ')[[1]][1])))})
	curr$lon <- sapply(pos,FUN=function(x){return(as.numeric(gsub('E','',strsplit(x,' ')[[1]][2])))})
	curr$datetime <- as.POSIXct(trimws(curr$timestamp), format='%d/%m/%Y %I:%M:%S %p')
	to.add <- data.frame(file=rep(file,nrow(curr)),timestamp=curr$timestamp,lat=curr$lat,lon=curr$lon,datetime=curr$datetime)
	
	#label with correct id, dyemark, name
	dyemarks <- sapply(meerkat.ids$dyemark,FUN=function(x){return(gsub('\\+','',x))})
	names <- as.character(meerkat.ids$name)
	ids <- as.character(meerkat.ids$id)
	filename <- strsplit(file,'/')[[1]] #filename (not full path)
	filename <- filename[length(filename)]
	parts <- strsplit(filename,'_')
	found <- F
	j <- 1
	while(!found){
		a <- c(grep(dyemarks[j],parts,ignore.case=T),grep(ids[j],parts,ignore.case=T),grep(names[j],parts,ignore.case=T))
		if(length(a)>0){
			found <- T
		} else{
			j <- j + 1
		}
	}
	to.add$device <- 'garmin'
	to.add$name <- names[j]
	to.add$dyemark <- dyemarks[j]
	to.add$id <- ids[j]
	
	#add to main data frame
	garmin.data <- rbind(garmin.data,to.add)
	
}

gps.data <- rbind(gps.data,garmin.data)

#clean up data
dates <- unique(as.Date(collar.data.good$datetime))
gps.data <- gps.data[which(as.Date(gps.data$datetime) %in% dates),]

#convert to utm
utm <- latlon.to.utm(LonsLat=cbind(gps.data$lon,gps.data$lat),southern_hemisphere=T,utm.zone='34')
gps.data$east <- utm[,1]
gps.data$north <- utm[,2]

#convert to local time by adding 2 hrs to gipsy5 data (garmin is already local)
gps.data$datetime[which(gps.data$device %in% c('collar','fluffy'))] <- gps.data$datetime[which(gps.data$device %in% c('collar','fluffy'))] + 2*60*60

#extract data only from the relevant dates and times
min.hours <- c(rep(9,15),rep(8,6))
max.hours <- c(rep(12,15),rep(11,6))
times <- data.frame()
day.start.idxs <- c(1)
for(i in 1:length(dates)){
	date <- dates[i]
	curr <- gps.data[which(as.Date(gps.data$datetime)==date),]
	hrs <- hour(curr$datetime)
	curr <- curr[which(hrs >= min.hours[i] & hrs < max.hours[i]),]
	times.curr <- data.frame(time=seq.POSIXt(from=min(curr$datetime),to=max(curr$datetime),by='sec'))
	times <- rbind(times,times.curr)
	day.start.idxs <- c(day.start.idxs,nrow(times)+1)
}

xs <- ys <- lats <- lons <- matrix(nrow=nrow(meerkat.ids),ncol=nrow(times))
for(i in 1:nrow(meerkat.ids)){
	curr <- gps.data[which(gps.data$id==as.character(meerkat.ids$id[i])),]
	idxs <- match(curr$datetime,times$time)
	non.na.idxs <- which(!is.na(idxs))
	xs[i,idxs[non.na.idxs]] <- curr$east[non.na.idxs]
	ys[i,idxs[non.na.idxs]] <- curr$north[non.na.idxs]
	lats[i,idxs[non.na.idxs]] <- curr$lat[non.na.idxs]
	lons[i,idxs[non.na.idxs]] <- curr$lon[non.na.idxs]
}
	
save(list=c('gps.data','xs','ys','lats','lons','times','day.start.idxs','meerkat.ids'),file='~/Dropbox/meerkats/data/Kalahari2017/RDATA/gps_level0.RData')




	
