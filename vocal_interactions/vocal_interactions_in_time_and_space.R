#This script uses cross-correlation functions to explore call/response dynamics as a function time and space

#It works as follows:
#Make a table of all calls within a given call category (set of calls), identify the one who gave the call as the 'caller'
#For each call, identify all other individuals ('responders') who were present (calls labeled) on that day, note their distance from the caller
#Also identify the distance between the caller and responder at the time when the call was given by the caller
#Save this table as callresp

#Then for each of row in this table, look at the call sequence of the 'responder' over time, lined up so that t = 0 = time of the 'caller' call
#Run a kernel over this time series to smooth it with a bandwidth bw = .1 sec
#Save the output to a matrix as callresp.seqs, such that the indices in the table match the indices in the matrix

#Finally, get the mean across all (or a subset) of rows of the matrix for different conditions:
# Different distance bins (aggregate only data from caller-responder pairs within a range of distances from one another)
# Different individuals + distance bins (aggregate only data from caller-responder pairs where the caller was a certain individual, and distance < 3 m since this was discovered to be the relevant range from the distance analysis)

#----------------PARAMETERS--------------------

#which year's data to use
year <- 2019 

#directory where data is stored
datadir <- '~/Dropbox/meerkats/meerkats_shared/data' 

#bandwidth of smoothing kernel (default 0.1)
bw <- .1 

#maximum time lag to consider (time since another individual called)
max.lag <- 30 

#time step to use for the sequence of times
step <- .01 

#list of call types to include in the set of calls by the initial caller (which determines the 0 point of the correlogram)
caller.calltypes <- c('cc','cchyb') 

#list of call types of include in the set of calls by the responder (determines the curve in the correlogram)
responder.calltypes <- c('cc','cchyb')

#distance bins to consider for the distance between caller and responder
dist.bins <- seq(0,12,3)

#whether to use self-responses instead of responses to others, for an analysis of individual calling periodicity (defaults to FALSE)
#normally this should be set to FALSE
self.responses <- FALSE 

#whether to use nonfocal calls within the same file as the 'caller' with the individual whose reocrding it is as 'respodner'
#this is mainly to check that the result is not driven by some weird synching issues, because we use the data frome the same recording
#normally this should be set to FALSE
use.nonfoc.as.triggers <- FALSE 


#---------------------SETUP-----------------------

#directories
audio.file <- paste(datadir, '/', year, '_ALL_CALLS_SYNCHED.csv',sep='')
gps.file <- paste(datadir ,'/', 'HM_COORDINATES_', year, '_sessions.RData', sep = '')
ind.file <- paste(datadir, '/', 'HM_INDIVIDUAL_INFO_', year,'.txt', sep = '')


#------------------LIBRARIES----------------------

#libraries
library(jcolors)
library(fitdistrplus)
library(viridis)
library(fields)

#------------------LOAD DATA-----------------------

#load audio data
load(gps.file)
calls.all <- read.csv(audio.file, header=T, sep='\t', stringsAsFactors=F)

#read in GPS data
load(gps.file)
ind.info <- read.csv(ind.file,sep='\t')

#--------------------PREPROCESS-------------------------------

#parse call types
#add a column for simple call type (parse calls in a sort of hierarchical way - should discuss this in more detail)
calls.all$callType <- tolower(calls.all$callType)
calls.all$callSimple <- 'oth'
calls.all$callSimple[which(calls.all$callType %in% c('s','sn','sx','sc','s?','snx'))] <- 's'
calls.all$callSimple[which(grepl('ag',calls.all$callType,ignore.case=T))] <- 'agg'
calls.all$callSimple[which(grepl('so',calls.all$callType,ignore.case=T))] <- 'soc'
calls.all$callSimple[which(grepl('chat',calls.all$callType,ignore.case=T))] <- 'chat'
calls.all$callSimple[which(grepl('mo',calls.all$callType,ignore.case=T))] <- 'mo'
calls.all$callSimple[which(grepl('ld',calls.all$callType,ignore.case=T))] <- 'ld'
calls.all$callSimple[which(grepl('lead',calls.all$callType,ignore.case=T))] <- 'ld'
calls.all$callSimple[which(grepl('al',calls.all$callType,ignore.case=T))] <- 'al'
calls.all$callSimple[which(grepl('cc',calls.all$callType,ignore.case=T))] <- 'cchyb'
calls.all$callSimple[which(calls.all$callType %in% c('cc','cc*','cc+','ccx','ccc','c'))] <- 'cc'
calls.all$callSimple[which(grepl('sync',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('beep',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('skip',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('start',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('stop',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('digging',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('eating',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('pause',calls.all$callType,ignore.case=T))] <- NA

#convert times to POSIXlt
calls.all$t0 <- as.POSIXlt(calls.all$t0GPS_UTC, tz = 'UTC')
calls.all$tf <- calls.all$t0 + calls.all$duration

#remove double-marked calls (same start and end time and same individual calling)
dups <- which(duplicated(cbind(calls.all$ind,as.numeric(calls.all$t0),as.numeric(calls.all$tf))))
calls.all <- calls.all[-dups,]

#get unique dates
dates <- unique(calls.all$date)


#---------------------------MAIN-------------------------------
#time sequence bins
tseq <- seq(-max.lag, max.lag, step)

#get together a table of all cc t0 and tf times from all individuals, as well as date and the start and end time of the labeled sequence
ccs.all <- data.frame()

for(d in 1:length(dates)){
  date <- dates[d]
  
  #get calls for that date
  calls.date <- calls.all[which(calls.all$date == date),]
  
  #get start and stop times for times when all are labeled - need to fix this a bit because sometimes there are multiiple start markers (e.g. VHMF010 on 20190712)
  starts <- calls.date[which(calls.date$callType %in% c('start')),c('ind','t0')]
  ends <- calls.date[which(calls.date$callType %in% c('stop','end')),c('ind','t0')]
  
  #for now just take the minimum 'start' for each individual as its start marker and the max end marker as its end time
  starts <- aggregate(starts$t0,by = list(starts$ind), min)
  ends <- aggregate(ends$t0,by = list(ends$ind), max)
  
  #find the latest start and earliest end times to be the group start and end time
  start.all <- max(starts$x)
  end.all <- min(ends$x)
  
  #make a simple table that just contains the cc's from each individual at that date and their t0 and tf times (as numeric)
  ccs <- calls.date[which(calls.date$callSimple %in% union(caller.calltypes, responder.calltypes) & !calls.date$nonFocal & !calls.date$unsureFocal & calls.date$isCall),]
  
  #if use.nonfoc.as.triggers is TRUE, do slightly differently
  if(use.nonfoc.as.triggers){
    ccs <- calls.date[which(calls.date$callSimple %in% union(caller.calltypes, responder.calltypes) & !calls.date$unsureFocal & calls.date$isCall),]
    ccs$ind[which(ccs$nonFocal==1)] <- paste(ccs$ind[which(ccs$nonFocal==1)],'NF',sep='_')
  }
  
  ccs <- ccs[,c('date','ind','callSimple','t0','tf')]
  ccs$t0 <- as.numeric(ccs$t0)
  ccs$tf <- as.numeric(ccs$tf)
  ccs$start <- as.numeric(start.all)
  ccs$end <- as.numeric(end.all)
  
  ccs.all <- rbind(ccs.all,ccs)
}

colnames(ccs.all) <- c('date','caller','callSimple','t0','tf','start','end')

#create a list of which individuals were labeled on which date
inds.present <- list()
for(i in 1:length(dates)){
  inds.present[[i]] <- unique(calls.all$ind[which(calls.all$date == dates[i])])
}

#create an expanded data frame that contains rows for each other 'responder' individual (not the caller)
callresp <- data.frame()
for(i in 1:nrow(ccs.all)){
  row <- ccs.all[i,]
  
  if(row$callSimple %in% caller.calltypes){
  
    date.idx <- which(dates == ccs.all$date[i])
    inds <- inds.present[[date.idx]]
    
    new.rows <- row[rep(1,length(inds)),]
    new.rows$responder <- inds
    
    callresp <- rbind(callresp, new.rows)
  }
  
}

#get distance between caller and responder at time of the call (t0)
callresp$distance <- NA
timeline.numeric <- as.numeric(timeLine, tz = 'UTC')
for(i in 1:nrow(callresp)){
  
  caller <- callresp$caller[i]
  responder <- callresp$responder[i]
  
  caller.idx <- match(caller,ind.info$code)
  responder.idx <- match(responder,ind.info$code)
  
  t0 <- round(callresp$t0[i]) ##
  t0.idx <- which(timeline.numeric==t0)
  
  if(length(t0.idx)!=0){
    callresp$distance[i] <- sqrt((allX[caller.idx,t0.idx] - allX[responder.idx,t0.idx])^2 + (allY[caller.idx,t0.idx] - allY[responder.idx,t0.idx])^2)
  }
  
}

#get rid of self-responses if self.responses==F, otherwise subset to only them
if(self.responses == T){
  callresp <- callresp[which(callresp$caller == callresp$responder),]
} else{
  callresp <- callresp[which(callresp$caller != callresp$responder),]
}

#create a matrix to hold call/response sequences
callresp.seqs <- matrix(data=NA,nrow=nrow(callresp),ncol = length(tseq))

#loop through and get individual call patterns starting at tf
for(i in 1:nrow(callresp.seqs)){
  zero.time <- callresp$t0[i]
  
  #if too close to end, don't use
  if((zero.time - max.lag > callresp$start[i]) & (zero.time + max.lag < callresp$end[i])){
    
    #get call times of the responder (only of the correct call type)
    call.times.foc <- ccs.all$t0[which(ccs.all$caller==callresp$responder[i] & (ccs.all$callSimple %in% responder.calltypes))]
    
    dt.foc <- call.times.foc - zero.time
    dt.foc <- dt.foc[which(abs(dt.foc) < max.lag)] #get rid of way early or way late calls to reduce compute time
    call.seq <- rep(0, length(tseq))
    
    if(length(dt.foc)>0){
      for(j in 1:length(dt.foc)){
        
        #for self-response version, don't include the current call
        if(self.responses & dt.foc[j] == 0){
          call.seq <- call.seq
        } else{
          call.seq <- call.seq + dnorm(tseq, mean = dt.foc[j], sd = bw)
        }
        
        
      }
    }
    
    callresp.seqs[i,] <- call.seq
    
  }
}

#---------------------PLOTTING-------------------

#ALL INDIVIDUALS - DIFFERENT DISTANCE RANGES

#Make a plot of cross-correlogram across all individuals

#collect the data
mean.call.rates <- matrix(NA,nrow=length(dist.bins)-1, ncol = length(tseq))
for(i in 2:length(dist.bins)){
  idxs <- which(callresp$distance >= dist.bins[i-1] & callresp$distance < dist.bins[i])
  mean.call.rates[i-1,] <- colMeans(callresp.seqs[idxs,],na.rm=T)
}

#if doing analysis on nonfocal calls within the same recording, do it differently
if(use.nonfoc.as.triggers){
  mean.call.rates <- rep(NA,length(tseq))
  idxs <- mapply(FUN=grepl,callresp$responder,callresp$caller)
  mean.call.rates <- colMeans(callresp.seqs[idxs,],na.rm=T)
  quartz(height = 8, width = 8)
  plot(tseq, mean.call.rates, lwd = 2,type='l',xlim=c(-10,10),xlab='Time lag (sec)', ylab = 'Call rate', cex.axis=1.5,cex.lab=1.5)
  abline(v=0,lty=2)
}

#make the plot
quartz(height = 12, width = 5)
par(mfrow=c(length(dist.bins)-1,1))
par(mar=c(6,5,1,1))
ymax <- max(mean.call.rates)+.01
#plot(NULL,xlim=c(-10,10),ylim=c(0,.15), xlab = 'Time lag (sec)', ylab = 'Call rate')
cols <- viridis(nrow(mean.call.rates))
for(i in 1:nrow(mean.call.rates)){
  plot(tseq, mean.call.rates[i,], col = cols[i], lwd = 2,type='l',xlim=c(-10,10),ylim=c(0,ymax),xlab='Time lag (sec)', ylab = 'Call rate', cex.axis=1.5,cex.lab=1.5)
  abline(v=0, lty=2)
  legend(x = 'topleft',legend=paste(dist.bins[i],'-', dist.bins[i+1],'m',sep=' '),col=cols[i],lwd=1.5, cex=1.5)
}


#SUBSET BY CALLER

#subset by caller (here only use 0-3 meter length scale as this is what was discovered before)
mean.call.rates <- matrix(NA,nrow=length(inds), ncol = length(tseq))
n.samps <- rep(NA, length(inds))
for(i in 1:length(inds)){
  idxs <- which(callresp$caller==inds[i] & callresp$distance < 3)
  mean.call.rates[i,] <- colMeans(callresp.seqs[idxs,],na.rm=T)
  n.samps[i] <- length(idxs)
}

quartz(height=12,width = 12)
par(mfrow=c(ceiling(length(inds)/3),3))
par(mar=c(6,5,1,1))

#set up color palette
cols <- jcolors('pal8')
ymax <- max(mean.call.rates,na.rm=T)+.01
for(i in 1:length(inds)){
  plot(tseq, mean.call.rates[i,], col = cols[i], lwd = 2,type='l',xlim=c(-5,5),ylim=c(0,ymax),xlab='Time lag (sec)', ylab = 'Call rate', cex.axis=1.5,cex.lab=1.5)
  abline(v=0, lty=2)
  legend(x = 'topleft',legend=inds[i],col=cols[i],lwd=1.5, cex=1.5)
}