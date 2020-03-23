#Extract call times for all individuals and files, and get them on (local) GPS time
#For now, fit a line to the synch calls and use this to interpolate (beeps only, for now)
#This is reasonable because drift does seem to be linear

#LIBRARIES
source('~/Dropbox/meerkats/code/general/meerkat_functions.R')
library(lubridate)
library(rhdf5)

#DIRECTORIES
call.dir <- '~/Dropbox/meerkats/data/Kalahari2017/labels_CSV'
synch.info.file <- '~/Dropbox/meerkats/data/Kalahari2017/RDATA/synch_info_all_2020-02-21.RData'

#LOAD MEERKAT IDs
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/ids.RData')

#LOAD TIMES
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/xy_level1.RData')

#FUNCTIONS
get.call <- function(instr){
	split.str <- strsplit(instr,split='/')[[1]]
	out <- split.str[1]
	return(as.character(out))
}
get.behav <- function(instr){
	split.str <- strsplit(instr,split='/')[[1]]
	out<-NA
	if(length(split.str)==2){
		out <- split.str[2]
	}
	return(as.character(out))
}

#MAIN

setwd(call.dir)

#Load synch information
load(synch.info.file)

#Get call data
calls <- data.frame()
files <- list.files()
for(i in 1:length(files)){ #for each file
  print(i)
	file <- files[i]
	print(file)
	
	#get the calls
	calls.curr <- read.table(file,header=T,sep='\t',stringsAsFactors=F,quote="",comment.char = "")
	calls.curr <- data.frame(call.name=calls.curr$Name,t0=calls.curr$Start,dur=calls.curr$Duration)
	
	#convert to seconds (rec time)
	calls.curr$t0.rec.sec <- sapply(calls.curr$t0,to.secs)
	calls.curr$tf.rec.sec <- sapply(calls.curr$dur,to.secs) + calls.curr$t0.rec.sec
	
	#convert to gps time
	
	#get all synch calls for that file - prioritize beeps if available
	synchs.curr <- synch.info[which(tolower(synch.info$filename)==tolower(file) & synch.info$synch.type=='beep' & !is.na(synch.info$speaker.time)),]
	print(synchs.curr)
	#if not enough beeps, go with synchs
	if(nrow(synchs.curr)<2){
	  synchs.curr <- synch.info[which(tolower(synch.info$filename)==tolower(file) & synch.info$synch.type=='synch' & !is.na(synch.info$speaker.time)),]
	}
	
	#fit a line
	x <- synchs.curr$rec.time
	y <- as.numeric(synchs.curr$gps.time)
	drift.mod <- lm(y~x)
	coefs <- drift.mod$coefficients
	out <- as.numeric(coefs[1]) + as.numeric(coefs[2])*calls.curr$t0.rec.sec
	calls.curr$t0.gps <- as.POSIXct(out,origin='1970-01-01 00:00.00 UTC',digits.secs=3)
	
	calls.curr$id <- synchs.curr$id[1]
	calls.curr$date <- synchs.curr$date[1]
	calls.curr$filename <- synchs.curr$filename[1]
	
	calls <- rbind(calls,calls.curr)
}

#Clean up call names

#create a new column to hold the call type
calls$call <- as.character(calls$call.name)

#replace default "Marker" with 'CC'
calls$call[grep('Marker',calls$call)] <- 'CC' 

#Synchs and beeps get replaced with NAs
calls$call[grep('SYNC',calls$call)] <- NA
calls$call[grep('BEEP',calls$call)] <- NA

#Remove trailing white spaces
calls$call <- trimws(calls$call)

#Separate calls, behaviors, nearest neighbor calls, markers
calls.only <- sapply(calls$call,FUN=get.call)
behavs.only <- sapply(calls$call,FUN=get.behav)
calls$call <- calls.only
calls$behavior <- behavs.only
calls$call[which(calls$call=='B')] <- NA
calls$nn <- NA
calls$nn[grep('NN',calls$call)] <- calls$call[grep('NN',calls$call)]
calls$call[grep('NN',calls$call)] <- NA
calls$mark <- NA
calls$mark[grep('START',calls$call)] <- calls$call[grep('START',calls$call)]
calls$mark[grep('END',calls$call)] <- calls$call[grep('END',calls$call)]
calls$mark[grep('RECON',calls$call)] <- calls$call[grep('RECON',calls$call)]
calls$mark[grep('RECOFF',calls$call)] <- calls$call[grep('RECOFF',calls$call)]
calls$mark[grep('INRANGE',calls$call)] <- calls$call[grep('INRANGE',calls$call)]
calls$mark[grep('OUTRANGE',calls$call)] <- calls$call[grep('OUTRANGE',calls$call)]
calls$mark[grep('PAUSE',calls$call)] <- calls$call[grep('PAUSE',calls$call)]
calls$mark[grep('RESUME',calls$call)] <- calls$call[grep('RESUME',calls$call)]
calls$mark[grep('[@]',calls$call)] <- calls$call[grep('[@]',calls$call)]

#mark chewing as a behavior
chews <- grep('CHEW',calls$call)
calls$behavior[chews] <- 'CHEW'
calls$call[chews] <- NA

calls$call[which(!is.na(calls$mark))] <- NA

#separate focal and nonfocal calls
calls$foc <- TRUE
nonfoc.calls <- grep('NON',calls$call)
unsure.calls <- grep('[*]',calls$call)
calls$foc[nonfoc.calls] <- FALSE
calls$foc[unsure.calls] <- NA
calls$foc[which(is.na(calls$call))] <- NA

#calls where you are not sure if it's a call or not
calls$unsure.call <- F
calls$unsure.call[grep('[#]',calls$call.name)] <- T

#calls where you are unsure about the type
calls$unsure.type <- F
calls$unsure.type[grep('[?]',calls$call.name)] <- T

#categorize as noisy or not
calls$noisy <- 0
calls$noisy[grep('X',calls$call)] <- 1
calls$call <- sub('X','',calls$call)

#categorize into basic categories
calls$call.type <- NA
lead.ccs <- grep('LEAD CC',calls$call)
ccs <- grep('CC',calls$call)
sns <- grep('SN',calls$call)
aggress <- unique(c(grep('AG',calls$call),grep('AGRESS',calls$call),grep('GROWL',calls$call)))
moves <- grep('MO',calls$call)
leads <- union(grep('LEAD',calls$call),grep('LD',calls$call))
alarms <- unique(c(grep('ALERT',calls$call),grep('ALARM',calls$call)))
socials <- grep('SOC',calls$call)
hybrids <- grep('HYB',calls$call)
hyb.ccs <- intersect(hybrids,ccs)
nas <- which(is.na(calls$call))
chats <- grep('CHAT',calls$call)
submits <- grep('SUBM',calls$call)
unks <- union(grep('UNK',calls$call),grep('OTH',calls$call))


calls$call.type <- 'OTH'
calls$call.type[unks] <- 'UNK'
calls$call.type[ccs] <- 'CC'
calls$call.type[sns] <- 'SN'
calls$call.type[moves] <- 'MOV'
calls$call.type[leads] <- 'LD'
calls$call.type[alarms] <- 'ALARM'
calls$call.type[aggress] <- 'AGG'
calls$call.type[socials] <- 'SOC'
calls$call.type[chats] <- 'CHAT'
calls$call.type[submits] <- 'SUBM'
calls$call.type[hybrids] <- 'HYB'
calls$call.type[hyb.ccs] <- 'HYB CC'
calls$call.type[lead.ccs] <- 'HYB CC'
calls$call.type[nas] <- NA

#Get the id, name, and dyemark of the focal individual
meerkat.names <- calls$id

#Get id and dye mark of focal
calls$id <- calls$dyemark <- calls$name <- calls$ind.idx <- NA
for(i in 1:nrow(meerkat.ids)){
	id <- as.character(meerkat.ids$id[i])
	dyemark <- as.character(meerkat.ids$dyemark[i])
	dyemark <- sub('\\+','',dyemark)
	name <- as.character(meerkat.ids$name[i])
	matches <- unique(c(grep(paste('_',id,'_',sep=''),calls$filename,ignore.case=T),grep(paste('_',dyemark,'_',sep=''),calls$filename,ignore.case=T),grep(paste('_',name,'_',sep=''),calls$filename,ignore.case=T)))
	if(length(matches)>0){
		calls$id[matches] <- id
		calls$dyemark[matches] <- dyemark
		calls$name[matches] <- name
		calls$ind.idx[matches] <- i
	}
}

#remove factors (replace with characters)
factors <- sapply(calls, is.factor)
calls[factors] <- lapply(calls[factors], as.character)

#-----------Add a column of time indexes to calls data frame
call.time.idxs <- match(as.character(calls$t0.gps),as.character(times$time)) #get indexes (of gps time) associated with calls
calls$gps.time.idx <- call.time.idxs #store these indexes in our calls table


save(list=c('calls'),file='~/Dropbox/meerkats/data/Kalahari2017/RDATA/calls_all_2020-02-21.RData')

#colnames(calls) <- gsub('\\.','_',colnames(calls))

h5createFile(file='~/Dropbox/meerkats/data/Kalahari2017/hdf5/calls_all_2020-02-21.h5')
h5write(file='~/Dropbox/meerkats/data/Kalahari2017/hdf5/calls_all_2020-02-21.h5',name='/calls',obj=calls)
H5close()




	