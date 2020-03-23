#ANIMATION WITH CALLS

#Out directory
base.dir <- '~/Dropbox/meerkats/animations/collective_movement'

#load level 0 gps data
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/gps_level0.RData')

#load calls data and colors
load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/calls_all.RData')
load('~/Dropbox/meerkats/data/Kalahari2017/RData/call_colors.R')

#Source functions
source('~/Dropbox/meerkats/code/general/meerkat_functions.R')

#number of meerkats and number of days
n.days <- length(day.start.idxs)-1
n.meerkats <- nrow(meerkat.ids)

calls$t.idx <- match(as.character(calls$t0.gps),as.character(times$time))

#make a plot of meerkat tracks with calls overlaid

t0 <- quantile(calls$t0.gps,.1)
tf <- quantile(calls$t0.gps,.97)

idx0 <- which(as.character(times$time)==as.character(t0))
idxf <- which(as.character(times$time)==as.character(tf))

#idx0 <- 156000
#idxf <- 156600

#set up for function
start.time <- idx0
end.time <- idxf
colors <- rainbow(n.meerkats)
tail.time <- 10
timesvec <- times$time
on.map <- F
plot.legend <- T
show.scale.bar <- F
ind.names <- meerkat.ids$id
show.scale.bar <- T
scale.bar.len <- 10
scale.bar.text.offset <- 0
utm.zone <- 36
zoom <- 16

ind.names <- c('Dom F (VLF206)','Dom M (VCVM001)','Sub F (VHMF001)','Sub M (VHMM002)','Sub M (VHMM003)','Sub M (VHMM006)','Sub M (VHMM007)')
trajectories.movie(base.dir=base.dir,inds=1:n.meerkats,plot.legend=T,start.time=start.time,end.time=end.time,step=1,lats=lats,lons=lons,on.map=on.map,tail.time=tail.time,times=timesvec,calls=calls,call.types=call.types,call.persist.time=10,ind.names=ind.names,show.scale.bar=show.scale.bar,scale.bar.len=scale.bar.len,scale.bar.text.offset=scale.bar.text.offset,utm.zone=utm.zone,zoom=zoom,scale.bar.text='10 m')

ts <- seq(idx0,idxf)
x <- xs[,ts]
y <- ys[,ts]

cols<-rainbow(n.meerkats)

xlims <- c(quantile(x,.001,na.rm=T),quantile(x,.999,na.rm=T))
ylims <- c(quantile(y,.001,na.rm=T),quantile(y,.999,na.rm=T))

for(i in 1:length(ts)){
	filename <- paste(outdir,'/',i,'.png',sep='')
	png(filename=filename,height=800,width=800)
	plot(NULL,xlim=xlims,ylim=ylims)
	points(xs[,ts[i]],ys[,ts[i]],pch=19,col=cols[i])
	dev.off()
}