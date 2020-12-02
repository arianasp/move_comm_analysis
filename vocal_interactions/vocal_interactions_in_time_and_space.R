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

#directory where code is stored
codedir <- '~/Dropbox/code_ari/move_comm_analysis/'

#bandwidth of smoothing kernel (default 0.1)
bw <- .2 

#maximum time lag to consider (time since another individual called)
max.lag <- 30 

#time step to use for the sequence of times
step <- .05 

#list of call types to include in the set of calls by the initial caller (which determines the 0 point of the correlogram)
caller.calltypes <- c('cc','cchyb') 

#list of call types of include in the set of calls by the responder (determines the curve in the correlogram)
responder.calltypes <- c('cc','cchyb')

#distance bins to consider for the distance between caller and responder
dist.bins <- seq(0,15,5)

#whether to use self-responses instead of responses to others, for an analysis of individual calling periodicity (defaults to FALSE)
#normally this should be set to FALSE
self.responses <- FALSE

#whether to use nonfocal calls within the same file as the 'caller' with the individual whose reocrding it is as 'respodner'
#this is mainly to check that the result is not driven by some weird synching issues, because we use the data frome the same recording
#normally this should be set to FALSE
use.nonfoc.as.triggers <- FALSE 

#whether to instead perform an analysis of whether an individual repeats its call in a certain time window after its initial call
#as a function of whether there have been any calls from nearby individuals (within repeat.dist.thresh meters) in the intervening time
#all sequences after a call are removed if there has been a call by anyone in the past repeat.time.thresh seconds
repeat.self.analysis <- TRUE
repeat.time.thresh <- 10
repeat.dist.thresh <- 10

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

#source functions
source(paste(codedir,'/general/meerkat_functions.R',sep=''))

#--------------------PREPROCESS-------------------------------

#add a column for simple call type (parse calls in a hierarchical way - see function description
calls.all$callSimple <- parse.meerkat.call.labels.to.simple.types(calls.all$callType)

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
  starts <- calls.date[which(calls.date$callType %in% c('start','START')),c('ind','t0')]
  ends <- calls.date[which(calls.date$callType %in% c('stop','end','STOP','END')),c('ind','t0')]
  
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

#use self responses if self.responses == T or if performing repeat self analysis (repeat.self.analysis==T), otherwise get rid of self responses
if(self.responses == T | repeat.self.analysis == T){
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
  if(((zero.time - max.lag) > callresp$start[i]) & ((zero.time + max.lag) < callresp$end[i])){
    
    #get call times of the responder (only of the correct call type)
    call.times.foc <- ccs.all$t0[which(ccs.all$caller==callresp$responder[i] & (ccs.all$callSimple %in% responder.calltypes))]
    
    dt.foc <- call.times.foc - zero.time
    dt.foc <- dt.foc[which(abs(dt.foc) < max.lag)] #get rid of way early or way late calls to reduce compute time
    call.seq <- rep(0, length(tseq))
    
    if(length(dt.foc)>0){
      for(j in 1:length(dt.foc)){
        
        #for self-response version, don't include the current call
        if((self.responses | repeat.self.analysis) & dt.foc[j] == 0){
          call.seq <- call.seq
        } else{
          call.seq <- call.seq + dnorm(tseq, mean = dt.foc[j], sd = bw)
        }
        
        
      }
    }
    
    callresp.seqs[i,] <- call.seq
    
  }
}

#if performing an analysis of whether you repeat a call depending on whether others have replied, need to compute some extra things
if(repeat.self.analysis == T){
  
  #----SETUP---
  #first the time of the previous call (by anyone within repeat.dist.thresh of the caller), the caller id, and the call type
  callresp$prevcall.t0 <- NA
  callresp$prevcall.caller <- NA
  callresp$prevcall.type <- NA
  callresp$n.inds.in.range <- NA
  
  #time of the next call (by anyone within repeat.dist.thresh of the caller), the caller id, and the call type
  callresp$nextcall.t0 <- NA
  callresp$nextcall.caller <- NA
  callresp$nextcall.type <- NA
  
  for(i in 1:nrow(callresp)){
    
    #get time of the focal call
    t0 <- callresp$t0[i] 
    
    #time index in timeLine
    t0.round <- round(t0)
    t.idx <- which(timeline.numeric == t0.round)
    
    #get individual index
    ind.idx <- which(ind.info$code == callresp$caller[i])
    
    #distances to all other individuals present
    dists <- sqrt((allX[,t.idx] - allX[ind.idx,t.idx])^2 + (allY[,t.idx] - allY[ind.idx,t.idx])^2)
    
    #distance to self is NA
    dists[ind.idx] <- NA
    
    #get the individuals that are in range (wihtin repeat.dist.thresh)
    inds.in.range.idxs <- which(dists < repeat.dist.thresh)
    
    #store number of individuals in range
    callresp$n.inds.in.range[i] <- length(inds.in.range.idxs)
    
    #if there are some indivudals in range, find the most recent call from any of them of any type, and next call of a type in responder.calltypes  
    if(length(inds.in.range.idxs)>0){
      
      #get individuals in range at the time of the focal call
      inds.in.range <- indInfo$code[inds.in.range.idxs]
      
      #get all calls from the individuals in range
      calls.inds.in.range <- calls.all[which(calls.all$ind %in% inds.in.range & calls.all$isCall==1 & calls.all$unsureFocal==0 & calls.all$nonFocal==0),]
      
      #get most recent previous call from the individuals in range (of any type)
      calls.before <- calls.inds.in.range$t0[which(calls.inds.in.range$t0 < t0)]
      if(length(calls.before) > 0){
        prev.call.time.inds.in.range <- max(calls.before)
        callresp$prevcall.t0[i] <- as.numeric(prev.call.time.inds.in.range, tz = 'UTC')
        prev.call.inds.in.range.idx <- which(calls.inds.in.range$t0 == prev.call.time.inds.in.range)[1]
        callresp$prevcall.caller[i] <- calls.inds.in.range$ind[prev.call.inds.in.range.idx]
        callresp$prevcall.type[i] <- calls.inds.in.range$callSimple[prev.call.inds.in.range.idx]
      }
      
      #get next call from individuals in range (of types in responder.calltypes)
      correct.calls.after <- calls.inds.in.range$t0[which(calls.inds.in.range$t0 > t0 & calls.inds.in.range$callSimple %in% responder.calltypes)]
      if(length(correct.calls.after)>0){
        correct.calls.after.numeric <- as.numeric(correct.calls.after, tz = 'UTC')
        next.call.time.inds.in.range <- min(correct.calls.after.numeric)
        callresp$nextcall.t0[i] <- as.numeric(next.call.time.inds.in.range)
        next.call.inds.in.range.idx <- which(calls.inds.in.range$t0 == next.call.time.inds.in.range)[1]
        callresp$nextcall.caller[i] <- calls.inds.in.range$ind[next.call.inds.in.range.idx]
        callresp$nextcall.type[i] <- calls.inds.in.range$callSimple[next.call.inds.in.range.idx]
      }
    }
  }
  
  callresp$dt.prev <- callresp$t0 - callresp$prevcall.t0
  callresp$dt.next <- callresp$nextcall.t0 - callresp$t0
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

#if performing analysis of whether calls are repeated when a reply is vs isn't given, make this plot
if(repeat.self.analysis){
  ##bins to measure individual call rate over (after a focal call), get response sequences only after t0 = 0
  bins <- tseq[which(tseq >= 0)]
  
  #get only positive time lags
  callresp.seqs.repeat <- callresp.seqs[,which(tseq >= 0)]
  
  #loop over bins to compute call rates when there has or hasn't been a response
  callrate.resp <- callrate.nonresp <- callrate.nonbrs <- n.resp <- n.nonresp <- n.nonbr <- rep(NA, length(bins)-1)
  for(i in 2:length(bins)){
    response.idxs <- which(callresp$dt.next < bins[i] & callresp$dt.prev > repeat.time.thresh)
    nonresponse.idxs <- which(callresp$dt.next > bins[i] & callresp$dt.prev > repeat.time.thresh)
    noneighbor.idxs <- which(callresp$n.inds.in.range == 0)
    
    callrate.resp[i] <- mean(callresp.seqs.repeat[response.idxs,i],na.rm=T)
    callrate.nonresp[i] <- mean(callresp.seqs.repeat[nonresponse.idxs,i],na.rm=T)
    callrate.nonbrs[i] <- mean(callresp.seqs.repeat[noneighbor.idxs,i],na.rm=T)
    n.resp[i] <- length(response.idxs)
    n.nonresp[i] <- length(nonresponse.idxs)
    n.nonbr[i] <- length(noneighbor.idxs)
  }
  
  quartz(height = 12, width = 8)
  par(mfrow=c(2,1))
  plot(bins, callrate.resp, type='l', col = 'red', lwd = 2, xlab = 'Time after call (sec)', ylab = c('Call rate'), cex.lab = 1.5, main = paste('Min silence = ', repeat.time.thresh, ' | Max dist = ', repeat.dist.thresh), ylim = c(0,0.2))
  lines(bins, callrate.nonresp, col = 'blue', lwd = 2, lty = 1)
  lines(bins, callrate.nonbrs, col = 'black', lwd = 2, lty = 1)
  legend('topright', legend = c('reply', 'no reply','no neighbor'), col = c('red','blue','black'), lty = c(1,1), cex = 1.5)
  plot(bins, n.resp, type = 'l', lwd = 2, xlab = 'Time (sec)', ylab = 'Number of events', col = 'red', cex = 2, cex.lab = 1.5, ylim = c(0, max(n.resp + n.nonresp + n.nonbr,na.rm=T)))
  lines(bins,n.nonresp, lwd = 2, col = 'blue')
  lines(bins,n.nonbr, lwd = 2, col = 'black')
}