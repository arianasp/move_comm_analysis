#Make some initial plots of the GPS data

#load level 0 gps data
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/gps_level0.RData')

#load calls data and colors
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/calls_all.RData')
load('~/Dropbox/meerkats/data/Kalahari2017/RData/call_colors.R')

#number of meerkats and number of days
n.days <- length(day.start.idxs)-1
n.meerkats <- nrow(meerkat.ids)

#PLOT 1
#make a plot of all data over all days (days colored the same)
cols <- rainbow(n.days)

plot(NULL,xlim=c(quantile(xs,.005,na.rm=T),quantile(xs,.995,na.rm=T)),ylim=c(quantile(ys,.005,na.rm=T),quantile(ys,.995,na.rm=T)),asp=1)
for(i in 1:n.days){
	ts <- day.start.idxs[i]:(day.start.idxs[i+1]-1)
	for(j in 1:n.meerkats){
		points(xs[j,ts],ys[j,ts],col=cols[i],pch=19,cex=0.01)
	}
}

#PLOT 2
#calls plotted on top of trajectories

#get time indexes of calls
calls$t.idx <- match(as.character(calls$t0.gps),as.character(times$time))

#make a plot of meerkat tracks with calls overlaid
t0 <- quantile(calls$t0.gps,.1)
tf <- quantile(calls$t0.gps,.9)

idx0 <- which(as.character(times$time)==as.character(t0))
idxf <- which(as.character(times$time)==as.character(tf))

idx0 <- 156000
idxf <- 156600

ts <- seq(idx0,idxf)
x <- xs[,ts]
y <- ys[,ts]

cols<-rainbow(n.meerkats)

plot(NULL,xlim=c(quantile(x,.001,na.rm=T),quantile(x,.999,na.rm=T)),ylim=c(quantile(y,.001,na.rm=T),quantile(y,.999,na.rm=T)),asp=1)
for(i in 1:n.meerkats){
	points(x[i,],y[i,],col=cols[i],pch=19,cex=0.1)
	ind.calls <- calls[which(calls$id==meerkat.ids$id[i] & is.na(calls$nonfoc) & calls$t.idx >= idx0 & calls$t.idx <= idxf),]
	xc <- xs[i,ind.calls$t.idx]
	yc <- ys[i,ind.calls$t.idx]
	colc <- call.types$col[match(ind.calls$call.type,call.types$type)]
	pchc <- call.types$sym[match(ind.calls$call.type,call.types$type)]
	points(xc,yc,col=as.character(colc),pch=pchc,cex=0.5)
}


#ANIMATION






