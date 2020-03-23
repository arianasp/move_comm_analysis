library(plotrix)

centr.x <- colMeans(xs,na.rm=T)
centr.y <- colMeans(ys,na.rm=T)
cols=rainbow(dim(xs)[1])

par(mfrow=c(1,2))
plot(centr.x[t:(t+dt)],centr.y[t:(t+dt)],type='l',asp=1)
draw.circle(centr.x[t],centr.y[t],10)
points(xs[,t],ys[,t],col=cols,pch=19)

#get relative xs and ys
plot(xs.rel[,t],ys.rel[,t],col=cols,pch=19,asp=1)
abline(h=0)
abline(v=0)
draw.circle(0,0,10)

#make relative positioning plot
load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/xy_level1.RData')
load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/ids.RData')
load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/positions_rel_to_centroid.RData')
cols=rainbow(dim(xs)[1])
n.inds <- dim(xs)[1]
cols = rep('black',n.inds)
cols[which(meerkat.ids$dom & meerkat.ids$sex=='F')] <- '#FF0000'
cols[which(meerkat.ids$dom & meerkat.ids$sex=='M')] <- '#0000FF'
cols[which(!meerkat.ids$dom & meerkat.ids$sex=='F')] <- '#FFAA66'
cols[which(!meerkat.ids$dom & meerkat.ids$se=='M')] <- '#66AAFF'
cols = tim.colors(n.inds)

t <- 1600
dt <- 600*6
trange <- t:(t+dt)
ts.axes <- c(t+200,t+2000,t+3000)

quartz(bg='black')
plot(NULL,xlim=c(min(xs[,trange],na.rm=T),max(xs[,trange],na.rm=T)),ylim=c(min(ys[,trange],na.rm=T),max(ys[,trange],na.rm=T)),bg='black',asp=1)
lines(c(min(xs[,trange],na.rm=T),min(xs[,trange],na.rm=T)+50),c(min(ys[,trange],na.rm=T),min(ys[,trange],na.rm=T)),col='white',lwd=2)
text(min(xs[,trange],na.rm=T)+25,min(ys[,trange],na.rm=T)+5,labels=c('50 m'),col='white')
for(i in 1:n.inds){
	lines(xs[i,trange],ys[i,trange],col=cols[i],lwd=1)
}
lines(centr.x[trange],centr.y[trange],col='white',lwd=2)
i <- 2100
R <- 20
for(i in ts.axes){
	arrows(centr.x[i]-R*cos(angs[i]),centr.y[i]-R*sin(angs[i]),centr.x[i]+R*cos(angs[i]),centr.y[i]+R*sin(angs[i]),col='white',lwd=3,length=0.1)
	arrows(centr.x[i]-R*cos(angs[i]-pi/2),centr.y[i]-R*sin(angs[i]-pi/2),centr.x[i]+R*cos(angs[i]-pi/2),centr.y[i]+R*sin(angs[i]-pi/2),col='white',lwd=3,length=0)
	points(centr.x[i],centr.y[i],pch=19,cex=1,col='white')
}