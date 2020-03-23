#Analyses of turning influence in meerkats (and some analyses of call rate vs spatial positioning)

#TODO: Incorporate breaks between (right now they are ignored, probably makes little difference but should be fixed for completeness)
#TODO: Deal with short gaps between files in some of the audio (e.g. dominant female) - only relevant for analyses using audio

#------------------------DIRECTORIES------------------------
codedir <- '~/Dropbox/meerkats/code'
plotdir <- '~/Dropbox/meerkats/results/turning_leadership/'

#------------------------PARAMETERS------------------------
subsamp <- 10 #subsample to once every subsamp sec to save processing time
R.heading <- 20 #radius to use for computing troop headings spatially
R.follower.dyad <- 5 #spatial scale for follower motion in dyad analysis
R.follower.centroid <- 5 #spatial scale for centroid following analysis
round <- 'HM2017' #whether to use data from 2017 or 2019, can be 'HM2017' or 'HM2019'
save_plots <- F
call.rate.window <- 60 #window (sec) over which to compute call rate

run_test <- F #whether to run a test (make plots to see if the measure is working correctly)
run_centroid_analysis <- F #whether to run analysis of turning influence on group centroid
run_dyad_analysis <- F #whether to run analysis of turning influence for each dyad
run_call_rate_analyses <- T #whether to run analysis of call rate vs relative positioning and group-level properties

#------------------------PRELIMINARIES------------------------
#set working directory
if(round == 'HM2017'){
  setwd('~/Dropbox/meerkats/data/Kalahari2017/RDATA/') #HM2017 data
}
if(round =='HM2019'){
  setwd('~/Dropbox/meerkats/meerkats_shared/05.02.2020/') #2019 data
}

#source some functions
source(paste(codedir,'/audio_gps_processing/spatially_discretized_headings.R',sep=''))
source(paste(codedir,'/general/meerkat_functions.R',sep=''))
source(paste(codedir, '/general/individual_and_group_level_properties_functions.R', sep = ''))
source('~/Dropbox/meerkats/meerkats_shared/code/call_interactions_functions.R')

#libraries
library(viridis)

#------------------------LOAD DATA------------------------
#load xy and ids data
if(round == 'HM2017'){
  load('xy_level1.RData') #HM2017 data
  load('ids.RData') #HM2017 data
  load('babysit_periods.RData')
  
  #format meerkat.ids same way as 2019
  meerkat.ids$Sex <- meerkat.ids$sex
  meerkat.ids$status <- 'Subordinate'
  meerkat.ids$status[which(meerkat.ids$dom==TRUE)] <- 'Dominant'
  meerkat.ids$shortCode <- c('F206','M01','F01','M02','M03','M06','M07')
}

if(round == 'HM2019'){
  load('HM_2019_ALL_COORDINATES.RData')
  xs <- allX
  ys <- allY
  times <- data.frame(time = timeLine)
  meerkat.ids <- read.table(file='HM_INDIVIDUAL_INFO.txt', sep = '\t', header=T)
}

#------------------------SETUP------------------------
#colors for meerkat ids
dom.m <- '#0000FF'
dom.f <- '#FF0000'
sub.m <- '#66CCFF'
sub.f <- '#FFCC66'
juv <- '#666666'
meerkat.ids$color <- NA
meerkat.ids$color[which(meerkat.ids$sex=='M')] <- sub.m
meerkat.ids$color[which(meerkat.ids$sex=='F')] <- sub.f
meerkat.ids$color[which(meerkat.ids$status=='Dominant' & meerkat.ids$sex=='F')] <- dom.f
meerkat.ids$color[which(meerkat.ids$status=='Dominant' & meerkat.ids$sex=='M')] <- dom.m
meerkat.ids$color[which(meerkat.ids$status=='Juvenile')] <- juv

#remove babysitters
if(round == 'HM2017'){
  for(i in 1:nrow(babysit.periods)){
    ind <- babysit.periods$id.idx[i]
    d <- babysit.periods$date.idx[i]
    t.idxs <- day.start.idxs[d]:(day.start.idxs[d+1]-1)
    xs[cbind(rep(ind,length(t.idxs)),t.idxs)] <- NA
    ys[cbind(rep(ind,length(t.idxs)),t.idxs)] <- NA
  }
}

#subsample to save time
xs.sub <- xs[,seq(1,ncol(xs),subsamp)]
ys.sub <- ys[,seq(1,ncol(ys),subsamp)]
times$idx <- 1:nrow(times)
times.sub <- times[seq(1,ncol(ys),subsamp),]

#get centroid
centr.x <- colMeans(xs.sub,na.rm=T)
centr.y <- colMeans(ys.sub,na.rm=T)

#get spatial headings (based on PREVIOUS times not future times!)
spat.heads <- spatial.headings(x = centr.x, y = centr.y, R = R.heading, backward = T)

#------------------------TEST------------------------

if(run_test){
  #pick two individuals and a start time
  leader <- 1
  follower <- 2
  t0 <- 1500
  dt <- 10
  
  #test to make sure it's working
  scores <- rep(NA,ncol(xs.sub))
  for(t0 in 1:ncol(xs.sub)){
    scores[t0] <- turning.leadership.score(xs.sub, ys.sub, spat.heads, leader = leader, follower = follower, t0 = t0, dt = 60, R = R)
  }
  
  idx <- 500
  idxs <- which(!is.na(scores))
  t0 <- idxs[idx]
  t.idxs <- t0:(t0+dt)
  xlim <- c(min(xs.sub[,t.idxs],na.rm=T), max(xs.sub[,t.idxs],na.rm=T))
  ylim <- c(min(ys.sub[,t.idxs],na.rm=T), max(ys.sub[,t.idxs],na.rm=T))
  plot(NULL,xlim=xlim,ylim=ylim,asp=1, main = as.character(scores[t0]))
  for(i in 1:nrow(xs)){
    if(i != leader & i != follower){
      lines(xs.sub[i, t.idxs], ys.sub[i, t.idxs],col='gray')
    }
  }
  points(xs.sub[leader, t.idxs],ys.sub[leader, t.idxs],col=viridis(dt+1),pch=19)
  lines(xs.sub[leader, t.idxs],ys.sub[leader, t.idxs],col='black',asp=1)
  points(xs.sub[follower, t.idxs],ys.sub[follower, t.idxs],col=viridis(dt+1),pch=19)
  lines(xs.sub[follower, t.idxs], ys.sub[follower, t.idxs], col = 'red')
  
  arrows(x0=centr.x[t0], y0=centr.y[t0], x1=centr.x[t0] + cos(spat.heads[t0])*10, y1=centr.y[t0] + sin(spat.heads[t0])*10, col='gray', lwd=3)
  arrows(x0=xs.sub[follower, t0], y0=ys.sub[follower,t0], x1=xs.sub[follower, t0 + dt], y1=ys.sub[follower, t0 + dt], col='red', lwd=3)
  lines(x=c(xs.sub[follower, t0], xs.sub[leader, t0]), y = c(y0=ys.sub[follower,t0], ys.sub[leader, t0]), col='black', lwd=3, lty = 2)
  print(scores[t0])
}

#------------------------CENTROID ANALYSIS------------------------

if(run_centroid_analysis){

  #get turning leadership scores for each individual on the group centroid
  scores.lead.centr <- matrix(NA,nrow=nrow(xs.sub),ncol=ncol(xs.sub))
  for(leader in 1:nrow(xs.sub)){
    timestamp()
    print(leader)
    for(t0 in 1:ncol(xs.sub)){
      scores.lead.centr[leader,t0] <- turning.leadership.score(xs = xs.sub, ys = ys.sub, spat.heads = spat.heads, leader = leader, follower = 'centroid',t0 = t0, dt = NULL, R = R.follower.centroid)
    }
  }

  #Plot individual scores relative to pulling centroid
  means.centr <- apply(scores.lead.centr,1,mean,na.rm=T)
  
  #set up color scale
  min.means <- min(means.centr[idxs.to.plot],na.rm=T)
  max.means <- max(means.centr[idxs.to.plot],na.rm=T)
  zlim.max <- max(abs(min.means),max(max.means))
  
  cols <- colorRampPalette(c('purple','white','green'))
  ncols <- 256
  cols <- cols(ncols)
  cols.boundaries <- seq(-zlim.max, zlim.max, length.out = ncols)
  
  height = 7
  width = 3
  quartz(height = height, width = width)
  idxs.to.plot <- 1:nrow(xs)
  par(mfrow = c(length(idxs.to.plot),3),mar=c(0.5,0.5,0.1,0.1))
  for(i in idxs.to.plot){
    plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',main='', bty='n')
    text(0.5,0.5,meerkat.ids$shortCode[i], adj=c(0.5,0.5), cex = 2, col = meerkat.ids$color[i])
    hist(scores.lead.centr[i,], breaks=seq(-1,1,.1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='black',freq=F, xlim = c(-1,1),ylim=c(0,2))
    rect(-1,0,1,2,border=NA, col=adjustcolor(cols[which.min(cols.boundaries < means.centr[i])],alpha = 0.6))
    plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',main='', bty='n')
    text(0.5,0.5,as.character(round(means.centr[i],digits=2)), adj=c(0.5,0.5), cex = 2, col = 'black')
  }
  
  if(save_plots){
    dev.copy2pdf(file=paste('~/Dropbox/meerkats/results/turning_leadership/',round,'_turning_leadership_on_centroid_grouphead',R.heading,'m_followerR', R.follower, 'm.pdf', sep=''))
    dev.off()
  }
}

#------------------------DYAD ANALYSIS------------------------

if(run_dyad_analysis){
  #Get scores for each pair
  scores <- array(NA, dim=c(nrow(xs.sub),nrow(xs.sub),ncol(xs.sub)))
  for(t0 in 1:ncol(xs.sub)){
    print(t0)
    for(leader in 1:nrow(xs.sub)){
      for(follower in 1:nrow(xs.sub)){
        if(leader != follower){
          scores[leader,follower,t0] <- turning.leadership.score(xs = xs.sub, ys = ys.sub, spat.heads = spat.heads, leader = leader, follower = follower, t0 = t0, dt = NULL, R = R.follower.dyad)
        }
      }
    }
  }
  
  #plot histograms (dyadic) and color by means
  means <- apply(scores,c(1,2),mean,na.rm=T)
  tots <- apply(scores,c(1,2),function(x){return(sum(!is.na(x)))})
  idxs.to.plot <- which(rowSums(tots) > 1000)
  
  #set up color scale
  min.means <- min(means[idxs.to.plot,idxs.to.plot],na.rm=T)
  max.means <- max(means[idxs.to.plot,idxs.to.plot],na.rm=T)
  zlim.max <- max(abs(min.means),max(max.means))
  
  cols <- colorRampPalette(c('purple','white','green'))
  ncols <- 256
  cols <- cols(ncols)
  cols.boundaries <- seq(-zlim.max, zlim.max, length.out = ncols)
  
  height = 7
  width = 7
  quartz(height = height, width = width)
  par(mar=c(2,2,2,2))
  par(mfrow=c(length(idxs.to.plot)+1,length(idxs.to.plot)+1),mar=c(0.5,0.5,0.1,0.1))
  #label top row
  plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',main='', bty='n')
  for(i in idxs.to.plot){
    plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',main='', bty='n')
    text(0.5,0.5,meerkat.ids$shortCode[i], adj=c(0.5,0.5), cex = 2, col = meerkat.ids$color[i])
  }
  
  for(i in idxs.to.plot){
    plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',main='', bty='n')
    text(0.5,0.5,meerkat.ids$shortCode[i], adj=c(0.5,0.5), cex = 2, col = meerkat.ids$color[i])
    for(j in idxs.to.plot){
      if(i==j){
        plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',main='', bty='n')
      }
      else{
        hist(scores[i,j,], breaks=seq(-1,1,.1),main='',xlab='',ylab='',xaxt='n',yaxt='n',col='black',freq=F, xlim = c(-1,1),ylim=c(0,2))
        rect(-1,0,1,2,border=NA, col=adjustcolor(cols[which.min(cols.boundaries < means[i,j])],alpha = 0.6))
      }
    }
  }
  
  if(save_plots){
    dev.copy2pdf(file=paste('~/Dropbox/meerkats/results/turning_leadership/',round,'_turning_leadership_by_dyad_grouphead',R.heading,'m_followerR', R.follower, 'm.pdf', sep=''))
    dev.off()
  }
}

#------------------------CALL RATE ANALYSES------------------------
if(run_call_rate_analyses){
  
  #LOAD DATA
  #load in call data
  if(round == 'HM2017'){
    load(file = 'calls_all_2020-02-21.RData')
  } 
  if(round == 'HM2019'){
    
    #preprocess calls to same format as 2017 data (rough version, only using categories we care about)
    calls <- read.csv(file = 'HM_2019_all_calls_synched.csv', sep ='\t')
    calls$call.type <- NA
    calls$call.type[which(calls$callType %in% c('cc','s+cc','cc+agg'))] <- 'CC'
    calls$call.type[which(calls$callType %in% c('mov','mo','ld','mo+ld','mov+ld','mov+s','s+mov','s+mo','s+ld'))] <- 'MOV'
    calls$call.type[which(calls$callType %in% c('al','al+al','al?','alarm','bark'))] <- 'ALARM'
    calls$call.type[which(calls$callType %in% c('s'))] <- 'SN'
    calls$foc <- !calls$unsureFocal & !calls$nonFocal & calls$isCall
    calls$gps.time.idx <- match(round(as.POSIXct(calls$t0GPS)), timeLine)
    calls$ind.idx <- match(calls$ind,meerkat.ids$code)
    calls$id <- calls$ind
    calls$mark <- NA
    calls$mark[which(calls$callType %in% c('stop','STOP','end','END'))] <- 'END'
    calls$mark[which(calls$callType %in% c('start','START'))] <- 'START'
  }
  
  #PRE-PROCESSING
  
  if(round == 'HM2017'){
    #dates that current have all inds labeled
    dates <- c('2017-08-06','2017-08-24','2017-08-25','2017-09-03','2017-09-05') 
  }
  if(round == 'HM2019'){
    dates <- c('20190713')
  }
  
  #get recording periods
  #NOTE: small gaps between files may reamin - TODO - FIX THIS
  #NOTE: stop/start times need to be fixed / clarified across all files! (for now, if not labeled we will assume start at first labeled call and end at last)
  rec.periods <- data.frame(period.idx = seq(1,length(dates)), date = dates, start.idx = NA, end.idx = NA, ids = NA, starts = NA, stops = NA)
  for(i in 1:nrow(rec.periods)){
    curr.calls <- calls[which(calls$date == rec.periods$date[i]),]
    ids <- unique(curr.calls$id)
    starts <- stops <- rep(NA, length(ids))
    for(j in 1:length(ids)){
      id <- ids[j]
      curr.calls.id <- curr.calls[which(curr.calls$id == id),]
      start.time <- curr.calls.id$gps.time.idx[which(curr.calls.id$mark %in% c('START','start'))]
      stop.time <- curr.calls.id$gps.time.idx[which(curr.calls.id$mark %in% c('END','end','STOP','stop'))]
      
      #TODO: make sure all files have labeled start and stop times!
      if(length(start.time)==0){
        start.time <- min(curr.calls.id$gps.time.idx,na.rm=T)
      }
      if(length(stop.time)==0){
        stop.time <- max(curr.calls.id$gps.time.idx,na.rm=T)
      }
      
      
      if(length(start.time)>0){
        starts[j] <- min(start.time)
      }
      if(length(stop.time)>0){
        stops[j] <- max(stop.time)
      }
    }
    latest.start <- max(starts, na.rm=T)
    earliest.end <- min(stops, na.rm=T)
    
    rec.periods$start.idx[i] <- latest.start
    rec.periods$end.idx[i] <- earliest.end 
    
    rec.periods$ids[i] <- list(ids)
    rec.periods$starts[i] <- list(starts)
    rec.periods$stops[i] <- list(stops)
  }
  
  #get time indexes to use
  t.idxs <- c()
  for(i in 1:nrow(rec.periods)){
    t.idxs <- c(t.idxs,rec.periods$start.idx[i]:rec.periods$end.idx[i])
  }
  
  #compute centroid position
  centr.x <- colMeans(xs,na.rm=T)
  centr.y <- colMeans(ys,na.rm=T)
  
  #get (previous) headings at specified time indexes (otherwise NA)
  spat.heads.past <- spatial.headings(centr.x, centr.y, R = R.heading, backward = T)
  
  #get (future) headings at specified time indexes (otherwise NA)
  spat.heads.future <- spatial.headings(centr.x, centr.y, R = R.heading, backward = F)
  
  #get call rates for each call type for each individual at the corresponding times
  cc.rates <- movld.rates <- sn.rates <- alarm.rates <- mov.rates <- ld.rates <- matrix(nrow=nrow(xs),ncol=ncol(xs))
  for(i in 1:nrow(xs)){
    print(i)
    calls.ind <- calls[which(calls$ind.idx == i & calls$foc==TRUE),]
    for(r in 1:nrow(rec.periods)){
      ts.curr <- seq(rec.periods$start.idx[r] + call.rate.window, rec.periods$end.idx[r], 1)
      for(j in 1:length(ts.curr)){
        t <- ts.curr[j]
        cc.rates[i,t] <- sum((calls.ind$call.type=='CC') & (calls.ind$foc==TRUE) & (calls.ind$gps.time.idx >= t - call.rate.window) & (calls.ind$gps.time.idx <= t),na.rm=T)
        mov.rates[i,t] <- sum((calls.ind$call.type %in% c('MOV')) & (calls.ind$foc==TRUE) & (calls.ind$gps.time.idx >= t - call.rate.window) & (calls.ind$gps.time.idx <= t),na.rm=T)
        ld.rates[i,t] <- sum((calls.ind$call.type %in% c('LD')) & (calls.ind$foc==TRUE) & (calls.ind$gps.time.idx >= t - call.rate.window) & (calls.ind$gps.time.idx <= t),na.rm=T)
        movld.rates[i,t] <- sum((calls.ind$call.type %in% c('MOV','LD')) & (calls.ind$foc==TRUE) & (calls.ind$gps.time.idx >= t - call.rate.window) & (calls.ind$gps.time.idx <= t),na.rm=T)
        alarm.rates[i,t] <- sum((calls.ind$call.type=='ALARM') & (calls.ind$foc==TRUE) & (calls.ind$gps.time.idx >= t - call.rate.window) & (calls.ind$gps.time.idx <= t),na.rm=T)
        sn.rates[i,t] <- sum((calls.ind$call.type=='SN') & (calls.ind$foc==TRUE) & (calls.ind$gps.time.idx >= t - call.rate.window) & (calls.ind$gps.time.idx <= t),na.rm=T)    
      }
    }
  }
  
  #get group turning angles
  dheads <- spat.heads.future - spat.heads.past
  dheads[which(dheads > pi)] <- dheads[which(dheads > pi)] - 2*pi
  dheads[which(dheads < -pi)] <- dheads[which(dheads < -pi)] + 2*pi
  
  #get positions relative to group heading and centroid
  relpos <- get.rel.positions(xs = xs, ys = ys, centr.x = centr.x , centr.y = centr.y, headings = spat.heads.past)
  xs.rel <- relpos$xs.rel
  ys.rel <- relpos$ys.rel
  
  #indexes to times when no alarms were given, no movld calls were given, or both (or just no nas)
  non.na.idxs <- which(!is.na(colSums(alarm.rates)))
  non.alarm.idxs <- which(colSums(alarm.rates)==0)
  non.movld.idxs <- which(colSums(movld.rates)==0)
  non.both.idxs <- intersect(non.alarm.idxs, non.movld.idxs)
  
  #---------------------PLOTS----------------
  
  #CALL RATES VS POSITION IN GROUP
  
  #get mean call rate as a function of position - front/back
  quartz(width=7, height = 7)
  idxs.to.use <- non.alarm.idxs
  xbins <- seq(-20,20,10)
  mids <- xbins[1:(length(xbins)-1)] + diff(xbins)/2
  plot(NULL, xlim=range(xbins), ylim = c(0,2), xlab = 'Distance front - back (m)', ylab = 'Mean lead call rate (calls / min)', cex.lab = 1.5, cex.axis=1.5)
  for(ind in 1:nrow(xs)){
    dat <- bin.plot(xvals = xs.rel[ind,idxs.to.use], yvals = ld.rates[ind,idxs.to.use], xbins = xbins,test.stat = 'mean', error.bar.quantiles = 'boot', show.plot = F, min.data.per.bin = 10)
    lines(mids, dat$vals, col = meerkat.ids$color[ind], lwd = 1)
    points(mids, dat$vals, col=meerkat.ids$color[ind],pch=19, cex = dat$tots / 3000)
  }
  
  #get mean call rate as a function of position - left/right
  quartz(width=7, height = 7)
  idxs.to.use <- non.alarm.idxs
  xbins <- seq(-20,20,5)
  mids <- xbins[1:(length(xbins)-1)] + diff(xbins)/2
  plot(NULL, xlim=range(xbins), ylim = c(0,5), xlab = 'Distance left - right (m)', ylab = 'Mean move call rate (calls / min)', cex.lab = 1.5, cex.axis=1.5)
  for(ind in 1:nrow(xs)){
    dat <- bin.plot(xvals = ys.rel[ind,idxs.to.use], yvals = mov.rates[ind,idxs.to.use], xbins = xbins,test.stat = 'mean', error.bar.quantiles = 'boot', show.plot = F, min.data.per.bin = 10)
    lines(mids, dat$vals, col = meerkat.ids$color[ind], lwd = 2)
    points(mids, dat$vals, col=meerkat.ids$color[ind],pch=19, cex = dat$tots / 1500)
  }
  
  #get mean move/lead call rate as a function of position - front/back
  quartz(width=7, height = 7)
  idxs.to.use <- non.alarm.idxs
  xbins <- seq(-20,20,5)
  mids <- xbins[1:(length(xbins)-1)] + diff(xbins)/2
  plot(NULL, xlim=range(xbins), ylim = c(0,5), xlab = 'Distance front - back (m)', ylab = 'Mean move/lead call rate (calls / min)', cex.lab = 1.5, cex.axis=1.5)
  for(ind in 1:nrow(xs)){
    dat <- bin.plot(xvals = xs.rel[ind,idxs.to.use], yvals = movld.rates[ind,idxs.to.use], xbins = xbins,test.stat = 'mean', error.bar.quantiles = 'boot', show.plot = F, min.data.per.bin = 10)
    lines(mids, dat$vals, col = meerkat.ids$color[ind], lwd = 1)
    points(mids, dat$vals, col=meerkat.ids$color[ind],pch=19, cex = dat$tots / 3000)
  }
  
  #get mean sn call rate as a function of position - front/back
  quartz(width=7, height = 7)
  idxs.to.use <- non.alarm.idxs
  xbins <- seq(-20,20,5)
  mids <- xbins[1:(length(xbins)-1)] + diff(xbins)/2
  plot(NULL, xlim=range(xbins), ylim = c(0,5), xlab = 'Distance left - right (m)', ylab = 'Mean short note call rate (calls / min)', cex.lab = 1.5, cex.axis=1.5)
  for(ind in 1:nrow(xs)){
    dat <- bin.plot(xvals = ys.rel[ind,idxs.to.use], yvals = sn.rates[ind,idxs.to.use], xbins = xbins,test.stat = 'mean', error.bar.quantiles = NULL, show.plot = F, min.data.per.bin = 10)
    lines(mids, dat$vals, col = meerkat.ids$color[ind], lwd = 1)
    points(mids, dat$vals, col=meerkat.ids$color[ind],pch=19, cex = dat$tots / 3000)
  }
  
  #get mean call rate as a function of position - left/right
  quartz(width=7, height = 7)
  idxs.to.use <- non.alarm.idxs
  xbins <- seq(-20,20,5)
  mids <- xbins[1:(length(xbins)-1)] + diff(xbins)/2
  plot(NULL, xlim=range(xbins), ylim = c(0,5), xlab = 'Distance left - right (m)', ylab = 'Mean move/lead call rate (calls / min)', cex.lab = 1.5, cex.axis=1.5)
  for(ind in 1:nrow(xs)){
    dat <- bin.plot(xvals = ys.rel[ind,idxs.to.use], yvals = movld.rates[ind,idxs.to.use], xbins = xbins,test.stat = 'mean', error.bar.quantiles = 'boot', show.plot = F, min.data.per.bin = 10)
    lines(mids, dat$vals, col = meerkat.ids$color[ind], lwd = 2)
    points(mids, dat$vals, col=meerkat.ids$color[ind],pch=19, cex = dat$tots / 3000)
  }
  
  
  #CALL RATES VS GROUP PROPERTIES
  
  time.diff <- 5*60 #time difference to use for computing speed (5 min)
  dx.centr <- c(centr.x[time.diff:(length(centr.x))] - centr.x[1:(length(centr.x)-time.diff+1)], rep(NA, time.diff-1))
  dy.centr <- c(centr.y[time.diff:(length(centr.x))] - centr.y[1:(length(centr.x)-time.diff+1)], rep(NA, time.diff-1))
  speed.future <- sqrt(dx.centr^2 + dy.centr^2) / time.diff * 60 # in m / min
  dx.centr <- c(rep(NA, time.diff-1), centr.x[time.diff:(length(centr.x))] - centr.x[1:(length(centr.x)-time.diff+1)])
  dy.centr <- c(rep(NA, time.diff-1), centr.y[time.diff:(length(centr.x))] - centr.y[1:(length(centr.x)-time.diff+1)])
  speed.past <- sqrt(dx.centr^2 + dy.centr^2) / time.diff * 60 # in m / min
  dspeed <- speed.future - speed.past
  
  #plot number of close calls across entire group vs mean speed
  quartz(width = 6, height = 6)
  bin.plot(xvals = colSums(cc.rates)[non.alarm.idxs], yvals = speed.future[non.alarm.idxs], xbins = c(0,10,20,30,40,50,70), error.bar.quantiles = 75, min.data.per.bin = 100, xlab='Group close call rate (calls / min)', ylab='Group speed (m / min)', plot.on.image = T, h = c(5,2))
  
  #plot number of move/lead calls across entire group vs mean speed
  quartz(width = 6, height = 6)
  bin.plot(xvals = movld.rates[1,][non.alarm.idxs], yvals = speed.future[non.alarm.idxs], xbins = c(0,1,2,3,4,5,7,10), error.bar.quantiles = 75, min.data.per.bin = 100, xlab='Group move/lead call rate (calls / min)', ylab='Group speed (m / min)', plot.on.image = T, h =c(1,2))
  
  #compare speed for different combos of DF & others move/lead call rates (heat map plot)
  ind <- 1
  idxs.to.use <- non.alarm.idxs
  xbins <- c(0,1,2,5,10)
  ybins <- c(0,1,2,5,10)
  out <- hist.2d(x = colSums(movld.rates[-ind,idxs.to.use]), y = movld.rates[ind,idxs.to.use], xbins = xbins, ybins = ybins, z = speed.future[idxs.to.use], xlab = 'Group calls', ylab = 'DF calls', output_freqs=T, output_plot = F)
  min.data.points <- 20
  out$histo[out$tots < min.data.points] <- NA
  quartz(width=6, height = 6)
  image.plot(z=out$histo, x = seq(1,length(xbins)-1), y = seq(1,length(ybins)-1), xlab = 'Group call rate', ylab = 'Dom Fem call rate', col= viridis(256), xaxt='n', yaxt='n', cex.lab = 1.5, cex.axis = 1.5)
  axis(1, at = seq(1,length(xbins)-1), labels = paste(xbins[1:(length(xbins)-1)], '-', as.character(xbins[2:length(xbins)]-1)))
  axis(2, at = seq(1,length(ybins)-1), labels = paste(ybins[1:(length(ybins)-1)], '-', as.character(ybins[2:length(ybins)]-1)))
  
  #simple plot of number of close callers vs speed (separate out DF calling vs. not)
  ind <- 4
  quartz(width=7, height = 6)
  xbins <- seq(0,8)
  nodf.idxs <- intersect(non.movld.idxs, which(cc.rates[ind,]==0))
  out.nodf <- bin.plot(xvals = colSums(cc.rates[,nodf.idxs]>0), yvals = speed.future[nodf.idxs], xbins = xbins, error.bar.quantiles = 'boot', test.stat='mean', show.plot=F)
  df.idxs <- intersect(non.movld.idxs, which(cc.rates[ind,]>0))
  out.df <- bin.plot(xvals = colSums(cc.rates[,df.idxs]>0), yvals = speed.future[df.idxs], xbins = xbins, error.bar.quantiles = 'boot', test.stat='mean', show.plot=F)
  plot(NULL,xlim=range(xbins), ylim=c(0,11),xlab='Number of close callers', ylab = 'Group speed (m / min)', cex.lab=1.5, cex.axis=1.5)
  mids <- out.df$xbins[1:(length(out.df$xbins)-1)]
  #arrows(mids, out.nodf$uppers, mids, out.nodf$lowers, angle=90, length = 0.1, col='black',lwd=2, code = 3)
  #arrows(mids, out.df$uppers, mids, out.df$lowers, angle=90, length = 0.1, col='black',lwd=2, code = 3)
  points(mids, out.df$vals, col='red',pch=19, cex = 2)
  points(mids, out.nodf$vals, col='black',pch=19, cex = 2)
  legend('topright',legend = c('Dom F calling','Dom F not calling'),pch=c(19,19), col = c('red','black'), cex = c(1.3,1.3))
  
  #simple plot of number of move callers vs speed (separate out DF calling vs. not)
  ind <- 1
  quartz(width=7, height = 6)
  xbins <- seq(0,5)
  nodf.idxs <- intersect(non.alarm.idxs, which(movld.rates[ind,]==0))
  nodf.idxs <- non.alarm.idxs
  out.nodf <- bin.plot(xvals = colSums(movld.rates[,nodf.idxs]>0), yvals = speed.future[nodf.idxs], xbins = xbins, error.bar.quantiles = 'boot', test.stat='mean', show.plot=F)
  df.idxs <- intersect(non.alarm.idxs, which(movld.rates[ind,]>0))
  out.df <- bin.plot(xvals = colSums(movld.rates[,df.idxs]>0), yvals = speed.future[df.idxs], xbins = xbins, error.bar.quantiles = 'boot', test.stat='mean', show.plot=F)
  plot(NULL,xlim=range(xbins), ylim=c(0,10),xlab='Number of move/lead callers', ylab = 'Group speed (m / min)', cex.lab=1.5, cex.axis=1.5)
  mids <- out.df$xbins[1:(length(out.df$xbins)-1)]
  #arrows(mids, out.nodf$uppers, mids, out.nodf$lowers, angle=90, length = 0.1, col='black',lwd=2, code = 3)
  #arrows(mids, out.df$uppers, mids, out.df$lowers, angle=90, length = 0.1, col='black',lwd=2, code = 3)
  #points(mids, out.df$vals, col='red',pch=19, cex = 2)
  points(mids, out.nodf$vals, col='gray',pch=19, cex = 2)
  #legend('topleft',legend = c('Dom F calling','Dom F not calling'),pch=c(19,19), col = c('red','black'), cex = c(1.3,1.3))
  
  
  #compare speed for different combos of DF & others close call rates (heat map plot)
  ind <- 2
  xbins <- c(0,11,21,31,41,51,61)
  ybins <- quantile(cc.rates[ind,idxs.to.use],c(0,.5,0.75,1))
  idxs.to.use <- non.both.idxs
  out <- hist.2d(x = colSums(cc.rates[-ind,idxs.to.use]), y = cc.rates[ind,idxs.to.use], xbins = xbins, ybins = ybins, z = speed.future[idxs.to.use], xlab = 'Group calls', ylab = 'DF calls', output_freqs=T, output_plot = F)
  min.data.points <- 10
  out$histo[out$tots < min.data.points] <- NA
  quartz(width=6, height = 6)
  image.plot(z=out$histo, x = seq(1,length(xbins)-1), y = seq(1,length(ybins)-1), xlab = 'Group call rate', ylab = 'Dom Fem call rate', col= viridis(256), xaxt='n', yaxt='n', cex.lab = 1.5, cex.axis = 1.5)
  axis(1, at = seq(1,length(xbins)-1), labels = paste(xbins[1:(length(xbins)-1)], '-', as.character(xbins[2:length(xbins)]-1)))
  axis(2, at = seq(1,length(ybins)-1), labels = paste(ybins[1:(length(ybins)-1)], '-', as.character(ybins[2:length(ybins)]-1)))
  
  #get call rate vs mean dyadic distance
  quartz()
  idxs.to.use <- non.alarm.idxs
  mean.dyad.dist <- get.mean.dyadic.dist(xs = xs, ys = ys)
  bin.plot(mean.dyad.dist, colSums(sn.rates[,idxs.to.use]), xbins = c(0,3,6,9,15,50), test.stat = 'mean', xlab = 'Mean dyadic distance (m)',ylab = 'CC rate')
  
  #Group spread (sd along each axis) vs call rate
  spreadx <- apply(xs.rel,2,sd,na.rm=T)
  spready <- apply(ys.rel,2,sd,na.rm=T)
  #qbins <- seq(0,1,.2)
  #xbins <- quantile(spreadx,qbins,na.rm=T)
  #ybins <- quantile(spready, qbins, na.rm=T)
  xbins <- ybins <- c(0,3,5,8,10,50)
  idxs.to.use <- non.both.idxs
  out <- hist.2d(x = spreadx, y = spready, xbins = xbins, ybins = ybins, z = colSums(cc.rates[,idxs.to.use],na.rm=T), xlab = 'Front - back spread (m)', ylab = 'Left - right spread (m)', output_freqs=T, output_plot = F)
  min.data.points <- 10
  out$histo[out$tots < min.data.points] <- NA
  quartz(width=6, height = 6)
  image.plot(z=out$histo, x = seq(1,length(xbins)-1), y = seq(1,length(ybins)-1), xlab = 'Front - back spread (m)', ylab = 'Left - right spread (m)', col= viridis(256), xaxt='n', yaxt='n', cex.lab = 1.5, cex.axis = 1.5)
  axis(1, at = seq(1,length(xbins)-1), labels = paste(xbins[1:(length(xbins)-1)], '-', as.character(xbins[2:length(xbins)]-1)))
  axis(2, at = seq(1,length(ybins)-1), labels = paste(ybins[1:(length(ybins)-1)], '-', as.character(ybins[2:length(ybins)]-1)))
  
  
  
  #CLOSE CALL BALANCE VS TURNING
  left <- cc.rates * (ys.rel > 0)
  right <- cc.rates * (ys.rel < 0)
  left.norm <- log((left+1) / (rowMeans(cc.rates, na.rm=T)+1))
  right.norm <- log((right+1) / (rowMeans(cc.rates, na.rm=T)+1))
  
  #get sum of left and right calls
  leftsum <- colSums(left, na.rm=T)
  rightsum <- colSums(right, na.rm=T)
  
  #make plot
  df.left.idxs <- which(ys.rel[1,] > 5)
  df.right.idxs <- which(ys.rel[1,] < 5)
  idxs.to.use <- intersect(non.both.idxs, df.left.idxs)
  idxs.to.use <- non.both.idxs
  quartz()
  xbins <- c(0,6,11,16,21,31,41,51,61)
  ybins <- c(0,6,11,16,21,31,41,51,61)
  out <- hist.2d(x = leftsum[idxs.to.use], y = rightsum[idxs.to.use], z = dheads[idxs.to.use]>0, xbins = xbins, ybins = ybins, zlim=c(0, 1), output_freqs = T, output_plot = F)
  cols <- colorRampPalette(c('purple','black','green'))
  min.data.points <- 10
  out$histo[out$tots < min.data.points] <- NA
  image.plot(z=out$histo, x = seq(1,length(xbins)-1), y = seq(1,length(ybins)-1), xlab = 'RHS call rate', ylab = 'LHS call rate', col=cols(256), xaxt='n', yaxt='n', cex.lab = 1.5, cex.axis = 1.5, zlim=c(0,1))
  axis(1, at = seq(1,length(xbins)-1), labels = paste(xbins[1:(length(xbins)-1)], '-', as.character(xbins[2:length(xbins)]-1)))
  axis(2, at = seq(1,length(ybins)-1), labels = paste(ybins[1:(length(ybins)-1)], '-', as.character(ybins[2:length(ybins)]-1)))
  
  
  #get probability turn left as a function of left call right minus right call rate
  idxs.to.use <- non.both.idxs
  dleftright <- leftsum - rightsum
  bin.plot(xvals = dleftright[non.both.idxs],yvals =dheads[non.both.idxs], xbins = c(-40,-20,-10,0,10,20,40))
  
  #TURNING INFLUENCE VS POSITION AND CALL RATE
  
  #compute turning influence score on the group centroid
  if(recalculate.turn.scores){
    subsamp.turn.score <- 30
    turn.score.centr <- matrix(NA,nrow=nrow(xs),ncol=ncol(xs))
    R.follower.centroid <- 5
    #for(j in 1:length(t.idxs)){
    for(j in seq(1,ncol(xs),subsamp.turn.score)){
      #if(j %% 100 == 0){
      #  print(round(j / length(t.idxs)*100))
      #  timestamp()
      #}
      print(round(j/ncol(xs)*100))
      for(i in 1:nrow(xs)){
        turn.score.centr[i,j] <- turning.leadership.score(xs = xs, ys = ys, spat.heads = spat.heads.past, leader = i, follower = 'centroid', dt = NULL, t0 = j, R = R.follower.centroid)
        #turn.score.centr[i,t.idxs[j]] <- turning.leadership.score(xs = xs, ys = ys, spat.heads = spat.heads.past, leader = i, follower = 'centroid', dt = NULL, t0 = t.idxs[j], R = R.follower.centroid)
      }
    }
    save(list=c('turn.score.centr'),file='~/Dropbox/meerkats/data/Kalahari2017/RDATA/HM2017_turning_influence_on_centroid_sub30_R5.RData')
  } else{
    load(file='~/Dropbox/meerkats/data/Kalahari2017/RDATA/HM2017_turning_influence_on_centroid_R20.RData')
  }
  
  #get normalized xs.rel and ys.rel (by group spread)
  spread.x <- apply(xs.rel, 2, sd, na.rm=T)
  spread.y <- apply(ys.rel, 2, sd, na.rm=T)
  spread.x.mat <- matrix(rep(spread.x, each=nrow(xs)), nrow=nrow(xs),ncol=ncol(xs))
  spread.y.mat <- matrix(rep(spread.y, each=nrow(xs)), nrow=nrow(xs),ncol=ncol(xs))
  xs.rel.norm <- xs.rel / spread.x.mat
  ys.rel.norm <- ys.rel / spread.y.mat
  
  #turning influence vs position along front-back axis
  idxs.to.use <- 1:ncol(xs)
  xbins <- seq(-2,2,.8)
  means.all <- uppers.all <- lowers.all <- matrix(NA,nrow=nrow(xs),ncol=length(xbins)-1)
  quartz(width = 10, height = 5)
  mids <- xbins[1:(length(xbins)-1)] + diff(xbins)/2
  par(mfrow=c(1,nrow(xs)),mar=c(5,5,0.1,0.1))
  for(ind in 1:nrow(xs)){
    if(ind==1){
      plot(NULL,xlim=range(xbins), ylim=c(0,1), ylab = 'P (turning influence > 0)', xlab = 'Left - right position in group (m)', cex.axis=1.5, cex.lab=1.5)
    } else{
      plot(NULL,xlim=range(xbins), ylim=c(0,1), ylab = '', xlab = 'Left - right position in group (m)', cex.axis=1.5, cex.lab=1.5)
    }
      out <- bin.plot(xvals = xs.rel.norm[ind,idxs.to.use], yvals = turn.score.centr[ind,idxs.to.use]>0, xbins = xbins, test.stat='mean',error.bar.quantiles = 'clopper',min.data.per.bin = 100, show.plot = F)
    non.nas <- which(!is.na(out$uppers))
    
    #estimate more realistic clopper-perason intervals based on downsampling (divide sample size by 20 because ACF goes to zero around 600 sec, and we already subsampled to once every 30 sec)
    tots <- round(out$tots/20)
    vals <- out$vals
    uppers <- lowers <- rep(NA, length(tots))
    for(c in 1:length(tots)){
      confi <- binom.test(n= tots[c], x = round(tots[c]*vals[c]))$conf.int
      uppers[c] <- confi[2]
      lowers[c] <- confi[1]
    }
    
    polygon(c(mids[non.nas],rev(mids[non.nas])),c(uppers[non.nas], rev(lowers[non.nas])), col = adjustcolor(meerkat.ids$color[ind],alpha.f=0.1),border=NA)
    lines(mids, out$vals, lwd=2, col= meerkat.ids$color[ind])
    points(mids, out$vals, pch=19, col = meerkat.ids$color[ind], cex = out$tot / 1000)
    abline(h=0.5,lty=2)
    abline(v=0, lty = 2)
    text(x=-2, y = 1, adj=0, labels =meerkat.ids$id[ind], col = meerkat.ids$color[ind], cex = 1.5)
  }
  
  #position and call rates vs turn score
  ind <- 1
  idxs.to.use <- non.both.idxs
  min.data <- 50
  xbins <- seq(-20,20,4)
  ybins <- c(0,1,100)
  out <- hist.2d(x = xs.rel[ind,idxs.to.use], y = cc.rates[ind,idxs.to.use], z=turn.score.centr[ind,idxs.to.use]>0, xbins = xbins, ybins = ybins, output_freqs=T)
  out$histo[out$tots < min.data] <- NA
  image.plot(out$histo, x = out$xbins, y = 1:length(out$ybins), zlim=c(0,1), col=cols(256), xlab = 'Front - back (m)', ylab = 'Close call rate (calls / min)')
  
  #turning influence on centroid vs close call rate
  quartz(width = 7, height = 7)
  idxs.to.use <- non.both.idxs
  qbins <- c(0,0.5,.75,1)
  xax <- 1:(length(qbins)-1)
  plot(NULL, xlim = range(xax), ylim = c(0,1), xaxt='n', xlab = 'Close call rate (calls / min)', ylab = 'Prob (turning influence > 0)')
  for(ind in 1:nrow(xs)){
    xbins <- quantile(cc.rates[ind, idxs.to.use],qbins,na.rm=T)
    out <- bin.plot(xvals = cc.rates[ind,idxs.to.use], turn.score.centr[ind,idxs.to.use]>0, xbins = xbins, test.stat = 'mean', error.bar.quantiles = 'boot', show.plot=F)
    polygon(c(xax, rev(xax)), c(out$uppers, rev(out$lowers)), col = adjustcolor(meerkat.ids$color[ind], alpha.f=0.2), border=NA)
    lines(seq(1,(length(xbins)-1)), out$vals, col = meerkat.ids$color[ind])
    points(seq(1,(length(xbins)-1)), out$vals, col = meerkat.ids$color[ind], pch = 19, cex = out$tots / 3000)
  }
  abline(h=0.5, lty=2)
  axis(1, at = xax, labels = c('Low (<50%)','Mid (50%-75%)','High (75-100%)'))
  
  
  #turning influence on centroid vs move/lead call rate
  quartz(width = 7, height = 7)
  idxs.to.use <- non.alarm.idxs
  xbins <- c(0,1,2,20)
  xax <- 1:(length(xbins)-1)
  plot(NULL, xlim = range(xax), ylim = c(0,1), xaxt='n', xlab = 'Close call rate (calls / min)', ylab = 'Prob (turning influence > 0)')
  for(ind in 1:nrow(xs)){
    out <- bin.plot(xvals = movld.rates[ind,idxs.to.use], turn.score.centr[ind,idxs.to.use]>0, xbins = xbins, test.stat = 'mean', error.bar.quantiles = 'boot', show.plot=F)
    polygon(c(xax, rev(xax)), c(out$uppers, rev(out$lowers)), col = adjustcolor(meerkat.ids$color[ind], alpha.f=0.2), border=NA)
    lines(xax, out$vals, col = meerkat.ids$color[ind])
    points(xax, out$vals, col = meerkat.ids$color[ind], pch = 19, cex = out$tots / 3000)
  }
  abline(h=0.5, lty=2)
  axis(1, at = xax, labels = c('0','1','>1'))
  
  #turning influence vs position along front-back axis and calling behavior
  idxs.to.use <- non.alarm.idxs
  ind <- 1
  xbins <- seq(-20,20,5)
  ybins <- c(0,1,4,20)
  min.data <- 100
  out <- hist.2d(x = xs.rel[ind,idxs.to.use], y = movld.rates[ind,idxs.to.use], z = turn.score.centr[ind,idxs.to.use]>0, xbins = xbins, ybins = ybins, output_freqs = T)
  out$histo[out$tots < min.data] <- NA
  image.plot(out$histo, x = 1:length(xbins), y = 1:length(ybins), col = cols(256), zlim=c(0,1))
  
  
  #TIME SERIES PLOTS
  quartz()
  par(mfrow=c(3,nrow(rec.periods)), mar=c(0,0,0,0))
  for(i in 1:nrow(rec.periods)){
    t <- rec.periods$start.idx[i]:rec.periods$end.idx[i]
    plot(t, colSums(movld.rates[-1,t],na.rm=T), type='l', col='black',lwd=2, ylim = c(0,40))
    lines(t,movld.rates[1,t], type='l', col='red',lwd=2)
  }
  for(i in 1:nrow(rec.periods)){
    t <- rec.periods$start.idx[i]:rec.periods$end.idx[i]
    plot(t, colSums(cc.rates[,t],na.rm=T), type='l', col='blue',lwd=2, ylim = c(0,80))
  }
  for(i in 1:nrow(rec.periods)){
    t <- rec.periods$start.idx[i]:rec.periods$end.idx[i]
    plot(t, speed.past[t], type='l', col='gray',lwd=2, ylim=c(0,15))
  }
  
  #CROSS-CORRELATIONS - call rates vs group speed
  #uses time shifts as null model but I don't think this null model really makes sense... how to assess significance of cross-correlation functions?
  max.lag <- 60*10
  ts1 <- colSums(movld.rates,na.rm=T)
  ts2 <- speed.past
  ccf.dat <- ccf(ts1, ts2, na.action = na.pass, lag.max = max.lag)
  plot(ccf.dat$lag, ccf.dat$acf,type='l')
  n.shifts <- 100
  acf.null <- matrix(nrow=length(ccf.dat$lag), ncol=n.shifts)
  for(i in 1:n.shifts){
    print(i)
    #shift one of the time series
    shift <- round(runif(1, min = -1, max = 1)*max.lag)
    if(shift>0){
      shift <- shift + max.lag
      ts2.shift <- ts2[shift:(length(ts2))]
      ts1.shift <- ts1[1:(length(ts1)-shift+1)]
    } else{
      shift <- shift - max.lag
      shift <- -shift
      ts1.shift <- ts1[shift:(length(ts1))]
      ts2.shift <- ts2[1:(length(ts2)-shift+1)]
    }
    tmp <- ccf(ts1.shift, ts2.shift, na.action = na.pass, lag.max = max.lag)
    acf.null[,i] <- tmp$acf
    
  }
  
  uppers <- apply(acf.null, 1, function(x){return(quantile(x,0.975,na.rm=T))})
  lowers <- apply(acf.null, 1, function(x){return(quantile(x,0.255,na.rm=T))})
  plot(ccf.dat$lag, ccf.dat$acf, type='l',lwd=2)
  lines(ccf.dat$lag, uppers)
  lines(ccf.dat$lag, lowers)
  
  
  #TOWARD-AWAY PLOTS
  load('toward_away_data_20200221.RData')
  twd.away.dat$call.rate.cc <- cc.rates[cbind(twd.away.dat$i, twd.away.dat$t)]
  twd.away.dat$call.rate.movld <- twd.away.dat$call.rate.mov <- NA
  twd.away.dat$call.rate.movld <- movld.rates[cbind(twd.away.dat$i, twd.away.dat$t)]
  twd.away.dat$alarm.present <- colSums(alarm.rates[,twd.away.dat$t])
  
  twd.away.dat$cc.rate.quantile <- NA
  for(i in 1:nrow(xs)){
    curr.idxs <- which(twd.away.dat$j == i)
    percentile <- ecdf(twd.away.dat$call.rate.cc[curr.idxs])
    twd.away.dat$cc.rate.quantile[curr.idxs] <- percentile(twd.away.dat$call.rate.cc[curr.idxs])
  }
  
  leaders <- c(1:7)
  qbins <- c(0.5,1)
  idxs.to.use <- which(twd.away.dat$j %in% c(leaders) & twd.away.dat$i %in% leaders & twd.away.dat$alarm.present==F)
  ybins <- c(0,0.5,0.75,1)
  xbins <- quantile(twd.away.dat$init.dist[idxs.to.use],seq(0,1,.2), na.rm=T)
  xbins <- c(0,100)
  out <- hist.2d(x = twd.away.dat$init.dist[idxs.to.use], y = twd.away.dat$cc,rate.quantile[idxs.to.use], z = cos(twd.away.dat$ang[idxs.to.use])>0, xbins = xbins, ybins = ybins, output_freqs = T, output_plot=F)
  
  out$histo[out$tots < 50] <- NA
  image.plot(out$histo, x = 1:length(out$xbins), y = 1:length(out$ybins), col=cols(256), xlab = 'Initial distance',ylab = 'Close call rate', zlim = c(0.1,0.9))
  
  #CORRELOGRAMS
  
  #move/lead
  max.lag <- 120
  njitters <- 100
  calls.movld <- calls[which(calls$call.type %in% c('MOV','LD') & calls$foc==TRUE),]
  means.jitter <- matrix(NA,nrow=njitters,,ncol=max.lag*2+1)
  for(k in 1:njitters){
    out <-kernel.call.rates(trigger.times = calls.movld$gps.time.idx[which(calls.movld$ind.idx==1)], call.times = calls.movld$gps.time.idx[which(calls.movld$ind.idx %in% c(2,3,4,5,6,7))], max.lag=max.lag, kernel.size = 5, start = 1, end = ncol(xs), hop.time = 1, jitter=T, jitter.time = max.lag*2)
    means.jitter[k,] <- colMeans(out)
  }
  out <-kernel.call.rates(trigger.times = calls.movld$gps.time.idx[which(calls.movld$ind.idx==1)], call.times = calls.movld$gps.time.idx[which(calls.movld$ind.idx %in% c(2,3,4,5,6,7))], max.lag=max.lag, kernel.size = 5, start = 1, end = ncol(xs), hop.time = 1, jitter=F)
  means.data <- colMeans(out)
  uppers <- apply(means.jitter, 2, function(x){return(quantile(x,0.975))})
  lowers <- apply(means.jitter, 2, function(x){return(quantile(x,0.025))})
  quartz(width=7,height=7)
  plot(seq(-max.lag,max.lag),colMeans(out), type= 'l', ylim=c(0,0.04), xlab = 'Time since dom F call (s)', ylab = 'Group call rate', cex.lab=1.5,cex.axis=1.5,lwd=2)
  polygon(c(seq(-max.lag,max.lag),seq(max.lag,-max.lag,-1)),c(uppers,rev(lowers)), col='#00000033', border=NA)
  abline(v=0,col='red')
  
  #cc
  max.lag <- 120
  njitters <- 100
  calls.cc <- calls[which(calls$call.type %in% c('CC') & calls$foc==TRUE),]
  means.jitter <- matrix(NA,nrow=njitters,,ncol=max.lag*2+1)
  for(k in 1:njitters){
    print(k)
    out <-kernel.call.rates(call.times = calls.cc$gps.time.idx[which(calls.cc$ind.idx==1)], trigger.times = calls.cc$gps.time.idx[which(calls.cc$ind.idx %in% c(2,3,4,5,6,7))], max.lag=max.lag, kernel.size = 5, start = 1, end = ncol(xs), hop.time = 1, jitter=T, jitter.time = max.lag*2)
    means.jitter[k,] <- colMeans(out)
  }
  out <-kernel.call.rates(call.times = calls.cc$gps.time.idx[which(calls.cc$ind.idx==1)], trigger.times = calls.cc$gps.time.idx[which(calls.cc$ind.idx %in% c(2,3,4,5,6,7))], max.lag=max.lag, kernel.size = 5, start = 1, end = ncol(xs), hop.time = 1, jitter=F)
  means.data <- colMeans(out)
  uppers <- apply(means.jitter, 2, function(x){return(quantile(x,0.975))})
  lowers <- apply(means.jitter, 2, function(x){return(quantile(x,0.025))})
  quartz(width=7,height=7)
  plot(seq(-max.lag,max.lag),colMeans(out), type= 'l', ylim=c(0.016,.022), xlab = 'Time since group member call (s)', ylab = 'Dom F call rate', cex.lab=1.5,cex.axis=1.5,lwd=2)
  polygon(c(seq(-max.lag,max.lag),seq(max.lag,-max.lag,-1)),c(uppers,rev(lowers)), col='#00000033', border=NA)
  abline(v=0,col='red')
  
  #cc vs movelead
  max.lag <- 120
  njitters <- 100
  calls.cc <- calls[which(calls$call.type %in% c('CC') & calls$foc==TRUE),]
  calls.movld <- calls[which(calls$call.type %in% c('MOV','LD') & calls$foc==TRUE),]
  means.jitter <- matrix(NA,nrow=njitters,,ncol=max.lag*2+1)
  for(k in 1:njitters){
    print(k)
    out <-kernel.call.rates(trigger.times = calls.movld$gps.time.idx[which(calls.movld$ind.idx %in% c(2,3,4,5,6,7))], call.times = calls.cc$gps.time.idx[which(calls.cc$ind.idx %in% c(1))], max.lag=max.lag, kernel.size = 5, start = 1, end = ncol(xs), hop.time = 1, jitter=T, jitter.time = max.lag*2)
    means.jitter[k,] <- colMeans(out)
  }
  out <-kernel.call.rates(trigger.times = calls.movld$gps.time.idx[which(calls.movld$ind.idx %in% c(2,3,4,5,6,7))], call.times = calls.cc$gps.time.idx[which(calls.cc$ind.idx %in% c(1))], max.lag=max.lag, kernel.size = 5, start = 1, end = ncol(xs), hop.time = 1, jitter=F)
  means.data <- colMeans(out)
  uppers <- apply(means.jitter, 2, function(x){return(quantile(x,0.975))})
  lowers <- apply(means.jitter, 2, function(x){return(quantile(x,0.025))})
  quartz(width=7,height=7)
  plot(seq(-max.lag,max.lag),means.data, type= 'l', ylim=c(0,.021), xlab = 'Time since group member move/lead call (s)', ylab = 'Dom F close call rate', cex.lab=1.5,cex.axis=1.5,lwd=2)
  polygon(c(seq(-max.lag,max.lag),seq(max.lag,-max.lag,-1)),c(uppers,rev(lowers)), col='#00000033', border=NA)
  abline(v=0,col='red')
  
  #DO CLOSE CALLS DURING MOVE CALL EVENTS SLOW DOWN THE GROUP?
  
  #first look to see whether giving move calls is associated with giving fewer close calls at an individual level
  ind <- 1
  out <- hist.2d(x = cc.rates[ind,non.alarm.idxs], y = movld.rates[ind,non.alarm.idxs], xbins = c(0,1,3,10,100), ybins= c(0,1,2,4,100))
  image.plot(out$histo, xlab = 'Close call rate', ylab = 'Move / lead call rate', col=viridis(256))
  
  #Make a heat map of future speed vs number of individual giving close calls and number giving move / lead calls
  quartz()
  min.data <- 100
  idxs.to.use <- intersect(which(movld.rates[1,]==0) ,non.alarm.idxs)
  n.moves <- colSums(movld.rates[,idxs.to.use]>0,na.rm=T)
  n.ccs <- colSums(cc.rates[,idxs.to.use]>0)
  z <- speed.future[idxs.to.use]
  xbins <- seq(0,max(n.moves)+1,na.rm=T)
  ybins <- seq(0,max(n.ccs)+1,na.rm=T)
  out <- hist.2d(x = n.moves, y = n.ccs, xbins = xbins, ybins = ybins, z = z, output_freqs = T, output_plot = F )
  out$histo[out$tots<100] <-NA
  image.plot(z=out$histo, x = xbins, y = ybins, col=viridis(256),xlab = 'Number of move callers',ylab = 'Number of close callers', xaxt = 'n', yaxt = 'n', cex.lab=1.5)
  axis(side=1, at = xbins[1:(length(xbins)-1)]+diff(xbins)/2, labels = as.character(xbins[1:(length(xbins)-1)]), cex.axis = 1.5)
  axis(side=2, at = ybins[1:(length(ybins)-1)]+diff(ybins)/2, labels = as.character(ybins[1:(length(ybins)-1)]), cex.axis = 1.5)
  
  
  #Make a heat map of group future speed vs. move call 
  ind <- 1
  movld.group <- colSums(movld.rates[-ind,non.alarm.idxs],na.rm=T)
  cc.df <- cc.rates[ind,non.alarm.idxs]
  z = speed.future[non.alarm.idxs]
  #roughly - low, mid, high - get cc rate bins based on individual quantiles
  xbins <- c(0,1,3,5,10,25)
  #qbins <- c(0,c(0.75,1))
  #ybins <- quantile(cc.df,qbins,na.rm=T)
  ybins <- c(0,1,100) #any close calls vs no close calls
  out <- hist.2d(x = movld.group, y = cc.df, xbins = xbins, ybins = ybins, z = z, output_freqs = T, output_plot = F )
  plot(x = seq(1,length(xbins)-1), out$histo[,1], xlab = 'Group move call rate',ylab = 'Mean group speed', pch=19, cex = 2, cex.lab=1.5, cex.axis=1.5,xaxt='n', ylim = c(0,8))
  points(x = seq(1,length(xbins)-1), out$histo[,2], pch = 19, cex = 2, col = 'red')
  axis(1, at = seq(1,length(xbins)-1), labels = c('0','1-2','3-4','5-9', '>10'), cex.axis = 1.5)
  
  
}


