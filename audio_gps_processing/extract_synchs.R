#extract synch calls and where they line up with seconds within a file
#output a data frame of this information for all files within a specific directory (indir)
#save it to an output directory

#LIBRARIES
library(stringr)
library(lubridate)
source('~/Dropbox/meerkats/code/general/meerkat_functions.R')

#DIRECTORIES
indir <- '~/Dropbox/meerkats/data/Kalahari2017/labels_CSV'
outdir <- '~/Dropbox/meerkats/data/Kalahari2017/RDATA'
gps.synch.file <- '~/Dropbox/meerkats/data/Kalahari2017/RDATA/synch_info_all.csv'

#MAIN

#set working directory to indir
setwd(indir)

#get list of files
files <- list.files()

#Read in info about each file and store in a data frame

#filename
file.info <- data.frame(filename=files)
file.info$filename <- as.character(file.info$filename)

#split into sound focals and collar data
foc.files <- grep('SOUNDFOC',file.info$filename)
collar.files <- seq(1,length(files))[-foc.files]

#individual id
file.info$id <- sapply(file.info$filename,FUN=function(x){return(strsplit(x,"_")[[1]][2])})

#dates
file.info$date <- NA
file.info$date[foc.files] <- sapply(file.info$filename[foc.files],function(x){return(strsplit(x,'_')[[1]][4])})
file.info$date[foc.files] <- sapply(file.info$date[foc.files],function(x){return(paste(substring(x,1,4),'-',substring(x,5,6),'-',substring(x,7,8),sep=''))})
file.info$date[collar.files] <- sapply(file.info$filename[collar.files],function(x){return(substring(x,which(strsplit(x,'')[[1]]=='(')+1,which(strsplit(x,'')[[1]]=='(')+10))})
file.info$date[collar.files] <- sapply(file.info$date[collar.files],function(x){return(gsub('_','-',x))})

#Load in GPS to audio synch table
gps.synch.dat <- read.table(gps.synch.file,header=T,sep=',')
gps.synch.dat$gps.datetime <- as.POSIXct(paste(gps.synch.dat$Date,gps.synch.dat$GPS.Time,sep=' '),format='%m/%d/%y %H:%M:%S')
gps.synch.dat$speaker.sec <- sapply(gps.synch.dat$Speaker.Time,to.secs)
gps.synch.dat$zero.datetime <- gps.synch.dat$gps.datetime - gps.synch.dat$speaker.sec
gps.synch.dat$min.datetime <- as.POSIXct(paste(gps.synch.dat$Date,gps.synch.dat$t0,sep=' '),format='%m/%d/%y %H:%M:%S')
gps.synch.dat$max.datetime <- as.POSIXct(paste(gps.synch.dat$Date,gps.synch.dat$tf,sep=' '),format='%m/%d/%y %H:%M:%S')


#Get synch info from files - both beeps (more accurate) and synchs (less accurate)
synch.info <- data.frame()
for(i in 1:nrow(file.info)){
	file <- file.info$filename[i]
	print(file)
	dat <- read.table(file,header=T,sep='\t',stringsAsFactors=F,comment.char = "",quote = )
	beeps <- dat[grep('BEEP',dat$Name),]
	synchs <- dat[grep('SYNCH',dat$Name),]
	if(nrow(beeps)>0){
	  beeps$Type <- 'beep'
	}
	if(nrow(synchs)>0){
	  synchs$Type <- 'synch'
	}
  print(beeps)
  print(synchs)
	dat <- rbind(beeps,synchs)
	dat$synchtime <- sapply(dat$Name,FUN=function(x){return(substring(x,which(strsplit(x,'')[[1]]==':')[1]-1,which(strsplit(x,'')[[1]]==':')[1]+5))})
	dat$synchtime.sec <- sapply(dat$synchtime,FUN=to.secs)
	dat$rectime.sec <- sapply(dat$Start,FUN=to.secs)
	gps.synch.dat.curr <- gps.synch.dat[which(as.character(date(gps.synch.dat$gps.datetime))==file.info$date[i]),]
	dat$gps.time <- gps.synch.dat.curr$zero.datetime[1] + dat$synchtime.sec
	diff.clock.idxs <- which(dat$gps.time > gps.synch.dat.curr$max.datetime[1])
	if(length(diff.clock.idxs)>0){
		for(j in 1:length(diff.clock.idxs)){
			idx <- diff.clock.idxs[j]
			dat$gps.time <- gps.synch.dat.curr$zero.datetime[2] + dat$synchtime.sec
		}
	}
	synch.info.curr <- data.frame(filename=rep(file,nrow(dat)),date=rep(file.info$date[i],nrow(dat)),id=rep(file.info$id[i],nrow(dat)),speaker.time=dat$synchtime.sec,rec.time=dat$rectime.sec,gps.time=dat$gps.time,synch.type=dat$Type)
	synch.info.curr$drift <- synch.info.curr$rec.time - synch.info.curr$rec.time[1] - (synch.info.curr$gps.time - synch.info.curr$gps.time[1])
	synch.info <- rbind(synch.info,synch.info.curr)
}

#resids <- c()
#for(i in 1:nrow(file.info)){
#	synch.info.curr <- synch.info[which(synch.info$filename==file.info$filename[i] & synch.info$synch.type=='beep'),]
#	x <- as.numeric(synch.info.curr$gps.time-synch.info.curr$gps.time[1])/60/60 #in hours
#	y <- as.numeric(synch.info.curr$drift)
#	plot(x,y)
#	print(i)
#	print(lm(y~x)$coefficients[2])
#	resids <- c(resids,lm(y~x)$residuals)
#}

save(list=c('synch.info'),file=paste(outdir,'/synch_info_all_2020-02-21.RData',sep=''))
