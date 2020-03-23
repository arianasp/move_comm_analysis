#Make a data frame containing
# t: time (once every 30 sec, say)
# i: focal individual
# j: nonfocal individaul (caller)
# xi0, yi0: initial position of focal
# xj0, yj0: initial position of caller
# xi0.rel, yi0.rel: initial position of focal relative to group axes
# xj0.rel, yj0.rel: initial position of nonfocal (caller) relative to group axes
# xif, yif: final position of focal after it has moved a distance R ( = 5 m)
# cc.rate: close call rate of the caller in the past t = 60 sec
# ld.rate: lead call rate of caller in past t = 60 sec
# mv.rate: move call rate of caller in past t = 60 sec

library(pracma)
library(fields)

#PARAMETERS

R <- 5 #radius to consider for subsequent movement direction (in m)
call.t.win <- 60 #time window over which to compute calls
dt <- 60

params <- list()
params$R <- R
params$call.t.win <- call.t.win
params$dt <- dt


#SOURCE
source('~/Dropbox/meerkats/code/general/meerkat_functions.R')
source('~/Dropbox/code_ari/baboons_code/step_selection/step_selection_funcs.R')

#LOAD DATA
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/positions_rel_to_centroid_prevhead.RData')
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/xy_level1.RData')
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/ids.RData')
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/babysit_periods.RData')
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/calls_all_2020-02-21.RData')

#PRE-PROCESS
#remove babysitters
for(i in 1:nrow(babysit.periods)){
  ind <- babysit.periods$id.idx[i]
  d <- babysit.periods$date.idx[i]
  t.idxs <- day.start.idxs[d]:(day.start.idxs[d+1]-1)
  xs[cbind(rep(ind,length(t.idxs)),t.idxs)] <- NA
  ys[cbind(rep(ind,length(t.idxs)),t.idxs)] <- NA
}

#GET NUMBER OF INDS AND TIMES
n.inds <- dim(xs)[1]
n.times <- dim(xs)[2]


#CREATE DATA FRAME
ts <- seq(1,n.times,dt)

dat <- data.frame()
for(i in 1:n.inds){
  for(j in 1:n.inds){
    if(!(i==j)){
      dat <- rbind(dat,data.frame(t=ts,i=rep(i,length(ts)),j=rep(j,length(ts))))
    }
  }
}

#add columns for initial positions
dat$xi0 <- xs[cbind(dat$i,dat$t)]
dat$yi0 <- ys[cbind(dat$i,dat$t)]
dat$xj0 <- xs[cbind(dat$j,dat$t)]
dat$yj0 <- ys[cbind(dat$j,dat$t)] 

#also relative to group center and heading
dat$xi0.rel <- xs.rel[cbind(dat$i,dat$t)]
dat$yi0.rel <- ys.rel[cbind(dat$i,dat$t)]
dat$xj0.rel <- xs.rel[cbind(dat$j,dat$t)]
dat$yj0.rel <- ys.rel[cbind(dat$j,dat$t)] 

#get maximum and minimum time index on a given day
dat$tmin <- sapply(dat$t,FUN=function(x){return(tail(day.start.idxs[which(day.start.idxs<=x)],n=1))})
dat$tmax <- sapply(dat$t,FUN=function(x){return(day.start.idxs[which(day.start.idxs>x)[1]]-1)})

dat$tf <- NA
for(r in 1:nrow(dat)){
  xcurr <- xs[dat$i[r],dat$t[r]:dat$tmax[r]]
  ycurr <- ys[dat$i[r],dat$t[r]:dat$tmax[r]]
  dx <- xcurr - dat$xi0[r]
  dy <- ycurr - dat$yi0[r]
  dists <- sqrt(dx^2+dy^2)
  tfs <- which(dists >= R)
  if(length(tfs)>0){
    dat$tf[r] <- tfs[1] + dat$t[r] - 1
  }
}

#get final position and angle of movement
dat$xif <- xs[cbind(dat$i,dat$tf)]
dat$yif <- ys[cbind(dat$i,dat$tf)]

#get angle of movement
v1x <- dat$xif-dat$xi0
v1y <- dat$yif-dat$yi0
v2x <- dat$xj0-dat$xi0
v2y <- dat$yj0-dat$yi0
dat$ang <- acos((v1x*v2x + v1y*v2y) / (sqrt((v1x^2+v1y^2)) * sqrt((v2x^2+v2y^2))))

#get call rate
files <- unique(calls$filename)
dat$call.rate.cc <- dat$call.rate.mov <- dat$call.rate.ld <- NA
for(f in 1:length(files)){
  print(f)
  calls.curr <- calls[which(calls$filename==files[f]),]
  start <- calls.curr$gps.time.idx[which(calls.curr$call.name=='START')]
  end <- calls.curr$gps.time.idx[which(calls.curr$call.name=='END')]
  caller <- calls.curr$ind.idx[1]
  idxs.dat <- which((dat$j==caller) & (dat$t >= (start+call.t.win)) & (dat$t <= end))
  
  times.cc <- calls.curr$gps.time.idx[which(calls.curr$call.type=='CC' & calls.curr$foc==T)]
  times.ld <- calls.curr$gps.time.idx[which(calls.curr$call.type=='LD'  & calls.curr$foc==T)]
  times.mov <- calls.curr$gps.time.idx[which(calls.curr$call.type=='MOV' & calls.curr$foc==T)]
  
  for(r in 1:length(idxs.dat)){
    tmin <- dat$t[idxs.dat[r]]-call.t.win
    tmax <- dat$t[idxs.dat[r]]
    dat$call.rate.cc[idxs.dat[r]] <- sum(times.cc >= tmin & times.cc < tmax,na.rm=T) / call.t.win * 60
    dat$call.rate.ld[idxs.dat[r]] <- sum(times.ld >= tmin & times.ld < tmax,na.rm=T) / call.t.win * 60
    dat$call.rate.mov[idxs.dat[r]] <- sum(times.mov >= tmin & times.mov < tmax,na.rm=T) / call.t.win * 60
  }
}

dat$init.dist <- sqrt((dat$xi0-dat$xj0)^2 + (dat$yi0-dat$yj0)^2)

#get whether there is a lead, move call in range of a certain time (from any meerkat) - 1 min forward or back
dat$ld.move.present <- dat$alarm.present <- NA
for(r in 1:nrow(dat)){
  print(r)
  dat$ld.move.present[r] <- sum((calls$call.type=='MOV' | calls$call.type=='LD') & abs(dat$t[r] - calls$gps.time.idx)<120,na.rm=T) > 0
  dat$alarm.present[r] <- sum((calls$call.type=='ALARM') & abs(dat$t[r] - calls$gps.time.idx)<60,na.rm=T) > 0
}

#SAVE DATA
twd.away.dat <- dat
save(list=c('twd.away.dat','params','xs','ys','xs.rel','ys.rel','calls'),file = '~/Dropbox/meerkats/data/Kalahari2017/RDATA/toward_away_data_20200221.RData')

