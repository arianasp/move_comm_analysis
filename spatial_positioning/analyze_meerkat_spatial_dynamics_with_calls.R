#ANALYZE MEERKAT SPATIAL DYNAMICS WITH RESPECT TO CALLS

library(pracma)

#LOAD DATA (generated from get_data_meerkat_spatial_dynamics_with_calls.R)
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/toward_away_data_20200221.RData')
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/xy_level1.RData')

#PARAMS
call.type <- 'CC'
speed.window <- 120 #window to use for computing speed

#SOURCE FUNCTIONS
source('~/Dropbox/meerkats/code/general/meerkat_functions.R')
source('~/Dropbox/code_ari/baboons_code/step_selection/step_selection_funcs.R')
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/ids.RData')

#COLORS
meerkat.ids$color <- c('red','blue','#FFAA66','#66AAFF','#66AAFF','#66AAFF','#66AAFF')
meerkat.ids$lty <- c(1,1,1,1,2,3,4)
meerkat.ids$pch <- c(19,19,19,19,18,17,15)

#ANALYSES TO PERFORM
call.rate.spatial.maps <- T
positioning.by.ind.rose.plots <- F

#LOCATIONS OF CALLS

#CALL LOCATION MAPS

#Get x and y position relative to group heading and normalized by sd along each axis
x.sd <- matrix(rep(apply(xs.rel,1,sd,na.rm=T),each=nrow(xs.rel)),nrow=nrow(xs.rel),ncol=ncol(xs.rel))
y.sd <- matrix(rep(apply(ys.rel,1,sd,na.rm=T),each=nrow(ys.rel)),nrow=nrow(ys.rel),ncol=ncol(ys.rel))
xs.rel.norm <- xs.rel / x.sd
ys.rel.norm <- ys.rel / y.sd

#get dyadic distance vs. time, as well as nn dist and identity
n.inds <- dim(xs)[1]
n.times <- dim(xs)[2]
dyad.dists <- array(NA,dim=c(n.inds,n.inds,n.times))
for(i in 1:(n.inds-1)){
  for(j in (i+1):n.inds){
    curr.dists <- sqrt((xs[i,] - xs[j,])^2 + (ys[i,] - ys[j,])^2)
      dyad.dists[i,j,] <- curr.dists
      dyad.dists[j,i,] <- curr.dists
  }
}

#mean dyadic dist
mean.dyad.dists <- apply(dyad.dists,3,mean,na.rm=T)

#nearest neighbor dist
nn.dists <- apply(dyad.dists,c(1,3),min,na.rm=T)
nn.dists[which(is.infinite(nn.dists))] <- NA

#nearest neighbor id
nn.ids <- apply(dyad.dists,c(1,3),which.min.random)



#position of centroid
x.centr <- colMeans(xs,na.rm=T)
y.centr <- colMeans(ys,na.rm=T)
n.tracked <- colSums(!is.na(xs))
x.centr[which(n.tracked<(n.inds-1))] <- NA
y.centr[which(n.tracked<(n.inds-1))] <- NA

#Get data where we have call info available
starts <- calls[which(calls$mark=='START'),]
ends <- calls[which(calls$mark=='END'),]
dates <- unique(calls$date)
starts.agg <- aggregate(starts$t0.gps,by=list(starts$ind.idx,starts$date),FUN=min)
ends.agg <- aggregate(ends$t0.gps,by=list(ends$ind.idx,ends$date),FUN=max)
start.times.by.day <- aggregate(starts.agg$x,by=list(starts.agg$Group.2),FUN=max)
end.times.by.day <- aggregate(ends.agg$x,by=list(ends.agg$Group.2),FUN=min)

#for now exclude 8-24 due to synch issues that need to be fixed
start.times.by.day <- start.times.by.day[-which(start.times.by.day$Group.1=='2017-08-24'),]
end.times.by.day <- end.times.by.day[-which(end.times.by.day$Group.1=='2017-08-24'),]

#Get daily time intervals when all individuals (with the group) are tracked and have audio
daily.time.intervals <- data.frame(date=start.times.by.day$Group.1,start=start.times.by.day$x,end=end.times.by.day$x)
daily.time.intervals$start.idx <- match(as.character(daily.time.intervals$start),as.character(times$time))
daily.time.intervals$end.idx <- match(as.character(daily.time.intervals$end),as.character(times$time))

#all time indexes
t.idxs <- c()
day.start.idxs.curr <- c(1)
row.idxs <- which(!is.na(daily.time.intervals$start.idx) & !is.na(daily.time.intervals$end.idx))
for(i in row.idxs){
  new.ts <- seq(daily.time.intervals$start.idx[i],daily.time.intervals$end.idx[i],60)
  t.idxs <- c(t.idxs,new.ts)
  day.start.idxs.curr <- c(day.start.idxs.curr,length(t.idxs)+1)
}


#get calls of that type
calls.curr <- calls[which(calls$call.type==call.type),]
calls.curr <- calls[which(calls$call.type=='LD' | calls$call.type=='MOV'),]

#get all positions at the relevant times
xs.all.norm <- xs.rel.norm[,t.idxs]
ys.all.norm <- ys.rel.norm[,t.idxs]

#get positions at the times of calls of the given type
xs.calls <- xs.rel.norm[cbind(calls.curr$ind.idx,calls.curr$gps.time.idx)]
ys.calls <- ys.rel.norm[cbind(calls.curr$ind.idx,calls.curr$gps.time.idx)]

#Make a histogram of meerkat positions when they make a certain call
bin.size <- 0.25
r.bins <- seq(-2,2,bin.size)
hist.2d(x=xs.calls,y=ys.calls,xbins=seq(-2,2,.2),ybins=seq(-2,2,.2))


#get call rate vs time
t.win <- 60
n.inds <- dim(xs)[1]

#get call rate data
call.rate.dat <- data.frame(t.idx=rep(t.idxs,n.inds),i=rep(1:n.inds,each=length(t.idxs)),t=times$time[t.idxs])
call.rate.dat$call.rate.cc <- NA
call.rate.dat$call.rate.ld <- NA
call.rate.dat$call.rate.mv <- NA
call.rate.dat$call.rate.sn <- NA
for(i in 1:n.inds){
  print(i)
  
  #close calls
  call.time.idxs <- calls$gps.time.idx[which(calls$call.type=='CC' & calls$ind.idx==i & calls$foc==T)]
  dat.idxs <- which(call.rate.dat$i == i)
  for(t in 1:length(t.idxs)){
    call.rate.dat$call.rate.cc[dat.idxs[t]] <- get.call.rate.fast(call.time.idxs,t.idxs[t],t.win)
  }
  
  #lead calls
  call.time.idxs <- calls$gps.time.idx[which(calls$call.type=='LD' & calls$ind.idx==i & calls$foc==T)]
  dat.idxs <- which(call.rate.dat$i == i)
  for(t in 1:length(t.idxs)){
    call.rate.dat$call.rate.ld[dat.idxs[t]] <- get.call.rate.fast(call.time.idxs,t.idxs[t],t.win)
  }
  
  #move calls
  call.time.idxs <- calls$gps.time.idx[which(calls$call.type=='MOV' & calls$ind.idx==i & calls$foc==T)]
  dat.idxs <- which(call.rate.dat$i == i)
  for(t in 1:length(t.idxs)){
    call.rate.dat$call.rate.mv[dat.idxs[t]] <- get.call.rate.fast(call.time.idxs,t.idxs[t],t.win)
  }
  
  #sn calls
  call.time.idxs <- calls$gps.time.idx[which(calls$call.type=='SN' & calls$ind.idx==i & calls$foc==T)]
  dat.idxs <- which(call.rate.dat$i == i)
  for(t in 1:length(t.idxs)){
    call.rate.dat$call.rate.sn[dat.idxs[t]] <- get.call.rate.fast(call.time.idxs,t.idxs[t],t.win)
  }
}

call.rate.dat$call.rate.cc.norm <- NA
for(i in 1:n.inds){
  idxs <- which(call.rate.dat$i==i)
  call.rate.dat$call.rate.cc.norm[idxs] <- log((call.rate.dat$call.rate.cc[idxs]+1) / (mean(call.rate.dat$call.rate.cc[idxs]+1,na.rm=T)))
}

#distances from front center side
call.rate.dat$dist.centr <- sqrt(xs.rel[cbind(call.rate.dat$i,call.rate.dat$t.idx)]^2 + ys.rel[cbind(call.rate.dat$i,call.rate.dat$t.idx)]^2)
call.rate.dat$dist.front <- ys.rel[cbind(call.rate.dat$i,call.rate.dat$t.idx)]
call.rate.dat$dist.side <- abs(xs.rel[cbind(call.rate.dat$i,call.rate.dat$t.idx)])

#normalized distances
call.rate.dat$dist.centr.norm <- sqrt(xs.rel.norm[cbind(call.rate.dat$i,call.rate.dat$t.idx)]^2 + ys.rel.norm[cbind(call.rate.dat$i,call.rate.dat$t.idx)]^2)
call.rate.dat$dist.front.norm <- ys.rel.norm[cbind(call.rate.dat$i,call.rate.dat$t.idx)]
call.rate.dat$dist.side.norm <- abs(xs.rel.norm[cbind(call.rate.dat$i,call.rate.dat$t.idx)])

#nearest neighbor distance and identity and group spread
call.rate.dat$mean.dyad.dist <- mean.dyad.dists[call.rate.dat$t.idx]
call.rate.dat$nn.dist <- nn.dists[cbind(call.rate.dat$i,call.rate.dat$t.idx)]
call.rate.dat$nn.id <- nn.ids[cbind(call.rate.dat$i,call.rate.dat$t.idx)]

call.rate.dat$n.tracked <- n.tracked[call.rate.dat$t.idx]

#get group speed before and after, and turning angle
dx2 <- (x.centr[call.rate.dat$t.idx+speed.window] - x.centr[call.rate.dat$t.idx])
dy2 <- (y.centr[call.rate.dat$t.idx+speed.window] - y.centr[call.rate.dat$t.idx])
dx1 <- (x.centr[call.rate.dat$t.idx] - x.centr[call.rate.dat$t.idx-speed.window])
dy1 <- (y.centr[call.rate.dat$t.idx] - y.centr[call.rate.dat$t.idx-speed.window])
dx0 <- (x.centr[call.rate.dat$t.idx-speed.window] - x.centr[call.rate.dat$t.idx-speed.window*2])
dy0 <- (y.centr[call.rate.dat$t.idx-speed.window] - y.centr[call.rate.dat$t.idx-speed.window*2])
call.rate.dat$group.speed.after <- sqrt(dx2^2 + dy2^2 ) / speed.window * 60
call.rate.dat$group.speed.before <- sqrt( dx1^2 + dy1^2 ) / speed.window * 60
call.rate.dat$group.speed.before.before <- sqrt( dx0^2 + dy0^2 ) / speed.window * 60
call.rate.dat$speed.change <- call.rate.dat$group.speed.after - call.rate.dat$group.speed.before
call.rate.dat$group.ta <- acos((dx1*dx2 + dy1*dy2) / (sqrt(dx1^2+dy1^2)*sqrt(dx2^2+dy2^2)))

#get group call rates for each call type
out <- aggregate(call.rate.dat$call.rate.cc,by=list(call.rate.dat$t.idx),FUN=sum)
call.rate.dat$group.call.rate.cc[match(call.rate.dat$t.idx,out$Group.1)] <- out$x
out <- aggregate(call.rate.dat$call.rate.mv,by=list(call.rate.dat$t.idx),FUN=sum)
call.rate.dat$group.call.rate.mv[match(call.rate.dat$t.idx,out$Group.1)] <- out$x
out <- aggregate(call.rate.dat$call.rate.ld,by=list(call.rate.dat$t.idx),FUN=sum)
call.rate.dat$group.call.rate.ld[match(call.rate.dat$t.idx,out$Group.1)] <- out$x
out <- aggregate(call.rate.dat$call.rate.sn,by=list(call.rate.dat$t.idx),FUN=sum)
call.rate.dat$group.call.rate.sn[match(call.rate.dat$t.idx,out$Group.1)] <- out$x

#get number of inds calling for each call type
out <- aggregate(call.rate.dat$call.rate.cc,by=list(call.rate.dat$t.idx),FUN=function(x){return(sum(x>0))})
call.rate.dat$n.callers.cc[match(call.rate.dat$t.idx,out$Group.1)] <- out$x
out <- aggregate(call.rate.dat$call.rate.mv,by=list(call.rate.dat$t.idx),FUN=function(x){return(sum(x>0))})
call.rate.dat$n.callers.mv[match(call.rate.dat$t.idx,out$Group.1)] <- out$x
out <- aggregate(call.rate.dat$call.rate.ld,by=list(call.rate.dat$t.idx),FUN=function(x){return(sum(x>0))})
call.rate.dat$n.callers.ld[match(call.rate.dat$t.idx,out$Group.1)] <- out$x
out <- aggregate(call.rate.dat$call.rate.sn,by=list(call.rate.dat$t.idx),FUN=function(x){return(sum(x>0))})
call.rate.dat$n.callers.sn[match(call.rate.dat$t.idx,out$Group.1)] <- out$x
out <- aggregate(call.rate.dat$call.rate.mv + call.rate.dat$call.rate.ld,by=list(call.rate.dat$t.idx),FUN=function(x){return(sum(x>0))})
call.rate.dat$n.callers.mvld[match(call.rate.dat$t.idx,out$Group.1)] <- out$x

#CALL RATE VS SPACE PLOTS

#bins
r.bins <- seq(0,2,.5)
theta.bins <- seq(-pi,pi,pi/2)
side.bins <- seq(0,2,.2)
front.bins <- seq(-2,2,.4)
nnd.bins <- seq(0,15,3)
spread.bins <- seq(0,25,5)

#data to use
curr <- call.rate.dat
vals <- curr$call.rate.sn

#call rate rose plot

#heatmap of location within group for each individual
heatmap.circ(xs=xs.rel.norm[cbind(curr$i,curr$t.idx)],ys=ys.rel.norm[cbind(curr$i,curr$t.idx)],theta.bins = seq(-pi,pi,pi/4),r.bins=r.bins,vals=vals,summary.stat='mean',cols=viridis(1024))

#Plot mean vs binned x value
x.bins <- seq(0,5,1)
vals.x <- call.rate.dat$n.callers.mvld
vals.y <- abs(call.rate.dat$group.speed.after)
lab.x <- 'Number of move or lead callers'
lab.y <- 'Mean speed after (m/min)'

means <- ns <- medians <- uppers <- lowers <- sds <- meanabs <- rep(NA,length(x.bins)-1)
for(i in 1:(length(x.bins)-1)){
  idxs.tmp <- which(vals.x >= x.bins[i] & vals.x < x.bins[i+1])
  means[i] <- mean(vals.y[idxs.tmp],na.rm=T)
  uppers[i] <- quantile(vals.y[idxs.tmp],0.75,na.rm=T)
  lowers[i] <- quantile(vals.y[idxs.tmp],0.25,na.rm=T)
  ns[i] <- sum(!is.na(vals.y[idxs.tmp]))
}
quartz()
par(mar=c(5.1, 5.1, 4.1, 2.1))
mids <- x.bins[1:(length(x.bins)-1)] 
y <- means
plot(mids,y,pch=19,xlab=lab.x,ylab=lab.y,cex=log(ns)*.5,cex.lab=2,cex.axis=2,ylim=c(min(lowers),max(uppers)))
arrows(mids,y,mids,uppers,angle=90,len=0.1,lwd=2)
arrows(mids,y,mids,lowers,angle=90,len=0.1,lwd=2)


#subset  by identity, plot on same plot
#get call rate as a function of distance to side and distance to front
x.bins <- c(0,1,20)
vals.x <- call.rate.dat$call.rate.ld
vals.y <- call.rate.dat$group.ta
lab.x <- 'Normalized distance toward front'
lab.y <- 'Mean call rate (calls / min)'


medians <- means <-  uppers <- lowers <- ns <- matrix(NA,nrow=n.inds,ncol=length(x.bins)-1)
for(j in 1:n.inds){
  individual.idxs <- which(call.rate.dat$i==j)
  curr <- call.rate.dat[individual.idxs,]
  vals.x.curr <- vals.x[individual.idxs]
  vals.y.curr <- vals.y[individual.idxs]
  for(i in 1:(length(x.bins)-1)){
    idxs <- which(vals.x.curr >= x.bins[i] & vals.x.curr < x.bins[i+1])
    medians[j,i] <- median(vals.y.curr[idxs],na.rm=T)
    means[j,i] <- mean(vals.y.curr[idxs],na.rm=T)
    uppers[j,i] <- quantile(vals.y.curr[idxs],0.75,na.rm=T)
    lowers[j,i] <- quantile(vals.y.curr[idxs],0.25,na.rm=T)
    ns[j,i] <- length(idxs)
  }
}

#plot
quartz(width=8,height=6)
mids <- x.bins[1:(length(x.bins)-1)] + diff(x.bins)/2
y <- means
xlim <- c(min(x.bins),max(x.bins))
xlim[2] <- xlim[2] + (xlim[2]-xlim[1])*.1
plot(NULL,xlim=xlim,ylim=c(min(y,na.rm=T),max(y,na.rm=T)),xlab=lab.x,ylab=lab.y)
for(i in 1:n.inds){
  lines(mids,y[i,],lwd=2,col=meerkat.ids$color[i],lty=meerkat.ids$lty[i])
  points(mids,y[i,],pch=meerkat.ids$pch[i],col=meerkat.ids$color[i])
}
legend('topright',col=meerkat.ids$color,legend=c('Dom F','Dom M','Sub F1','Sub M2','Sub M3','Sub M6','Sub M7'),lty=meerkat.ids$lty,pch=meerkat.ids$pch)

#Remove data when lead and move calls were present
dat <- twd.away.dat
dat.normal <- dat[which(!is.na(dat$call.rate.cc) & !dat$alarm.present),]

#Plot call rate vs distance apart
quartz()

dist.bins <- quantile(dat.normal$init.dist,seq(0,1,length.out=6),na.rm=T)
means <- uppers<- lowers <- ns <- matrix(NA,nrow=n.inds,ncol=length(dist.bins)-1)
for(j in 1:n.inds){
  curr <- dat.normal[which(dat.normal$j==j),]
  for(i in 1:(length(dist.bins)-1)){
    idxs <- which(curr$init.dist >= dist.bins[i] & curr$init.dist < dist.bins[i+1])
    means[j,i] <- mean(curr$call.rate.cc[idxs],na.rm=T)
    uppers[j,i] <- sd(curr$call.rate.cc[idxs],na.rm=T)
    ns[j,i] <- sum(!is.na(curr$call.rate.cc[idxs]))
  }
}
quartz()
mids <- dist.bins[2:length(dist.bins)]-diff(dist.bins)/2
cols <- rainbow(n.inds)
par(mfrow=c(2,4))
for(i in 1:n.inds){
  plot(mids,means[i,],lwd=2,col=cols[i],type='l',xlab='Distance (m)',ylab='Call rate (calls/min)',cex.lab=1.5,main=meerkat.ids$name[i])
  points(mids,means[i,],pch=19,,col=cols[i])
}

#plot toward away as a function of call rate - first raw
call.rate.bins <- c(0,1,100)
dist.bins <- quantile(dat$init.dist,seq(0,1,length.out=6),na.rm=T)
quartz()
idxs <- which(dat$j>=0)
test <- hist.2d(x=dat.normal$call.rate.cc[idxs],y=dat.normal$init.dist[idxs],xbins=call.rate.bins,ybins=dist.bins,vals=dat.normal$ang[idxs]<(pi/2),output_freqs = T)
image.plot(test$heatmap,xlab='call rate',ylab='distance')

means <- ns <- rep(NA,(length(call.rate.bins)-1))
for(i in 1:(length(call.rate.bins)-1)){
  idxs <- which(dat$call.rate.mov >= call.rate.bins[i] & dat$call.rate.mov < call.rate.bins[i+1])
  means[i] <- mean(dat$ang[idxs]<(pi/4),na.rm=T)
  ns[i] <- sum(!is.na(dat$ang[idxs]))
}

for(i in 1:n.inds){
  idxs <- which(dat$j==i)
  all.rates <- dat$call.rate.cc[idxs]
  quant.bins <- quantile(all.rates,seq(0,1,.1),na.rm=T)
  dat$call.rate.cc.norm[idxs] <- (all.rates - min(all.rates,na.rm=T)) / max(all.rates,na.rm=T)
}

#POSITIONING BY INDIVIDUAL - ROSE PLOTS
if(positioning.by.ind.rose.plots){
  max.r <- 1.5
  n.rs <- 10
  r.bins.equalized <- c(0,sqrt(seq(1,n.rs)))
  r.bins.equalized <- r.bins.equalized / max(r.bins.equalized) * max.r
  theta.bins <- seq(-pi,pi,pi/8)
  zlim=c(0,.039)
  
  idx.starts <- day.start.idxs[c(1,11,16)]
  idx.ends <- day.start.idxs[c(11,16,22)]-1
  
  quartz(width=8,height=8)
  par(mfrow=c(3,3),mar=c(1,1,1,1))
  for(i in c(1,6,7)){
    for(d in 1:length(idx.starts)){
      heatmap.circ(xs=xs.rel.norm[i,idx.starts[d]:idx.ends[d]],ys=ys.rel.norm[i,idx.starts[d]:idx.ends[d]],theta.bins = theta.bins,r.bins=r.bins.equalized,vals=NULL,summary.stat='mean',cols=viridis(1024),zlim=zlim)
    }
  }
  
}

#TOWARD AWAY VS DISTANCE BY CALL RATE
callers <- c(1, 2, 3, 4, 5, 6, 7)
curr <- twd.away.dat[which(twd.away.dat$j %in% callers & !is.na(twd.away.dat$call.rate.mov) & twd.away.dat$alarm.present==F),]
#curr <- twd.away.dat[which(twd.away.dat$j %in% callers),]
quantiles <- seq(0,1,1/5)
dist.bins <- quantile(curr$init.dist,quantiles,na.rm=T)
callrate.bins <- quantile(curr$call.rate.mov + curr$call.rate.ld,na.rm=T)
callrate.bins <- c(0,1,100)
out <-hist.2d(x=curr$init.dist,y=curr$call.rate.mov,vals = curr$ang<=(pi/2),xbins = dist.bins,ybins = callrate.bins,output_freqs=T)

x <- dist.bins
quartz()
plot(NULL,xlim=c(min(x),30),ylim=c(40,100),xlab='Initial distance between meerkats (m)',ylab='% approach')

arrows(x[1:(length(x)-1)],out$heatmap[,1]*100,x[2:length(x)],out$heatmap[,1]*100,length=0,lwd=3,col='black')
arrows(x[1:(length(x)-1)],out$heatmap[,2]*100,x[2:length(x)],out$heatmap[,2]*100,length=0,lwd=3,col='red')
arrows(x[1:(length(x)-1)],out$heatmap[,3]*100,x[2:length(x)],out$heatmap[,3]*100,length=0,lwd=3,col='blue')

abline(h=50,lty=2,lwd=2)


#TOWARD VS AWAY ROSE PLOTS - MOVE AND LEAD CALLS
quartz()
callers <- c(2:7)
curr <- twd.away.dat[which(twd.away.dat$j %in% callers & (twd.away.dat$ld.move.present==F) & twd.away.dat$alarm.present==F),]
curr <- curr[which(curr$call.rate.cc>4),]
dx.rel <- curr$xj0.rel - curr$xi0.rel
dy.rel <- curr$yj0.rel - curr$yi0.rel
theta.rel <- atan2(dy.rel,dx.rel)
heatmap.circ(rs=curr$init.dist,thetas=theta.rel,r.bins=seq(0,20,5),theta.bins = seq(-pi,pi,pi/4),vals=curr$ang<pi/2,zlim = c(0,1))

