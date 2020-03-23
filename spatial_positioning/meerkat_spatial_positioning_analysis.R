#Look for consistent individual differences in where meerkats position themselves in the group
#(across days)

#PARAMS
min.n.tracked <- 4 #minimum number of individuals tracked to include data
R.spat <- 10 #radius to use for spatial discretization to determine troop heading
recompute.rel.pos <- T
make.hists <- F
prev.head <- T #whether the compute headings from previous direction (T) or from future direction (F)

#LIBRARIES
library(fields)
library(ICC)

#SOURCE FUNCTIONS
source('~/Dropbox/code_ari/baboons_code/step_selection/spatial_discretization.R')
source('~/Dropbox/meerkats/code/general/meerkat_functions.R')
source('~/Dropbox/meerkats/code/general/individual_and_group_level_properties_functions.R')

#LOAD DATA

if(recompute.rel.pos){
	load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/xy_level1.RData')
	load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/ids.RData')
	load('~/Dropbox/meerkats/data/Kalahari2017/RDATA/babysit_periods.RData')

	#REMOVE THE BABYSITTERS
	for(i in 1:nrow(babysit.periods)){
		ind <- babysit.periods$id.idx[i]
		d <- babysit.periods$date.idx[i]
		t.idxs <- day.start.idxs[d]:(day.start.idxs[d+1]-1)
		xs[cbind(rep(ind,length(t.idxs)),t.idxs)] <- NA
		ys[cbind(rep(ind,length(t.idxs)),t.idxs)] <- NA
	}
	

	#get troop centroid and number tracked
	xs.centr <- t(colMeans(xs,na.rm=T))
	ys.centr <- t(colMeans(ys,na.rm=T))
	n.tracked <- t(colSums(!is.na(xs)))

	#remove centroid data when not enough tracked
	xs.centr[which(n.tracked < min.n.tracked)] <- NA
	ys.centr[which(n.tracked < min.n.tracked)] <- NA

	#get headings by projecting R.spat m back and R.spat m into the future and drawing a line
	heads.x <- array(NA,dim=dim(xs.centr))
	heads.y <- array(NA,dim=dim(xs.centr))
	for(d in 1:(length(day.start.idxs)-1)){
		curr.x <- xs.centr[day.start.idxs[d]:(day.start.idxs[d+1]-1)]
		curr.y <- ys.centr[day.start.idxs[d]:(day.start.idxs[d+1]-1)]
		print(d)
		if(!prev.head){
			for(i in 1:length(curr.x)){
				x <- curr.x[i]
				y <- curr.y[i]
				idx.fwd <- 1
				found <- F
				#go forward to find distance R m away in future
				while(!found){
					if((i+idx.fwd) > length(curr.x)){
						x.after <- NA
						y.after <- NA
						found <- T
					} else{
						dist <- sqrt( (x-curr.x[i+idx.fwd])^2 + (y-curr.y[i+idx.fwd])^2)
						if(is.na(dist)){
							x.after <- NA
							y.after <- NA
							found <- T
						} else{
							if(dist > R.spat){
								x.after <- curr.x[i+idx.fwd]
								y.after <- curr.y[i+idx.fwd]
								found <- T
							}
						}
					}
					idx.fwd <- idx.fwd + 1
				}
				dx <- x.after - x
				dy <- y.after - y
				heads.x[day.start.idxs[d]+i-1] <- dx / sqrt(dx^2 + dy^2)
				heads.y[day.start.idxs[d]+i-1] <- dy / sqrt(dx^2 + dy^2)
			}
		} else{
			#computations for using previous data to compute heading
			for(i in 1:length(curr.x)){
				x <- curr.x[i]
				y <- curr.y[i]
				idx.bwd <- 1
				found <- F
				#go forward to find distance 5 m away in future
				while(!found){
					if((i-idx.bwd) < 1){
						x.before <- NA
						y.before <- NA
						found <- T
					} else{
						dist <- sqrt( (x-curr.x[i-idx.bwd])^2 + (y-curr.y[i-idx.bwd])^2)
						if(is.na(dist)){
							x.before <- NA
							y.before <- NA
							found <- T
						} else{
							if(dist > R.spat){
								x.before <- curr.x[i-idx.bwd]
								y.before <- curr.y[i-idx.bwd]
								found <- T
							}
						}
					}
					idx.bwd <- idx.bwd + 1
				}
				dx <- x - x.before
				dy <- y - y.before
				heads.x[day.start.idxs[d]+i-1] <- dx / sqrt(dx^2 + dy^2)
				heads.y[day.start.idxs[d]+i-1] <- dy / sqrt(dx^2 + dy^2)
			}
		}
	}

	#project individual positions into group reference frame

	#first get distance along x and y axis (world reference frame) from centroid
	xs.centr.mat <- matrix(rep(xs.centr,each=dim(xs)[1]),nrow=dim(xs)[1],ncol=dim(xs)[2])
	ys.centr.mat <- matrix(rep(ys.centr,each=dim(xs)[1]),nrow=dim(xs)[1],ncol=dim(xs)[2])
	dx <- xs - xs.centr.mat
	dy <- ys - ys.centr.mat

	#then rotate into group reference frame
	angs <- atan2(heads.y,heads.x) #angles
	rot.angs <- -angs + pi /2
	rot.angs.mat <- matrix(rep(rot.angs,each=dim(xs)[1]),nrow=dim(xs)[1],ncol=dim(xs)[2])
	xs.rel <- dx*cos(rot.angs.mat) - dy*sin(rot.angs.mat)
	ys.rel <- dx*sin(rot.angs.mat) + dy*cos(rot.angs.mat)
			
	#save relative positions and headings

	save(file='/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/positions_rel_to_centroid_prevhead.RData',list=c('xs.rel','ys.rel','angs','R.spat'))
} else{

load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/xy_level1.RData')
load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/positions_rel_to_centroid.RData')
load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/ids.RData')
}

#'independent' day ranges
day.ranges <- list()
day.ranges[[1]] <- 1:10
day.ranges[[2]] <- 11:15
day.ranges[[3]] <- 16:21

#Compute histograms
if(make.hists){
	png('/Users/astrandb/Dropbox/meerkats/results/spatial_positioning/positioning_hists2.png')
	par(mfrow=c(7,3),mar=c(0,0,0,0))
	for(j in 1:length(day.ranges)){
		n.inds <- dim(xs.rel)[1]
		for(i in 1:n.inds){
			x <- y <- c()
			for(k in day.ranges[[j]]){
				x <- c(x,xs.rel[i,day.start.idxs[k]:(day.start.idxs[k+1]-1)])
				y <- c(y,ys.rel[i,day.start.idxs[k]:(day.start.idxs[k+1]-1)])
				print(k)
			}
			hist.2d(x=x,y=y,xbins=seq(-20,20,1),ybins=seq(-20,20,1),axes=F,zlim=c(0,.02),imageplot=F)
			abline(h=0,lty=2,col='white')
			abline(v=0,lty=2,col='white')
		}
	}
	dev.off()
}

#Compute medians
n.inds <- dim(xs.rel)[1]
n.days <- length(day.start.idxs)-1
frac.front <- median.dist.centr <- median.abs.x <- median.abs.y <- matrix(NA,nrow=n.inds,ncol=n.days)
for(i in 1:n.inds){
	for(j in 1:n.days){
		xs.curr <- xs.rel[i,day.start.idxs[j]:(day.start.idxs[j+1]-1)]
		ys.curr <- ys.rel[i,day.start.idxs[j]:(day.start.idxs[j+1]-1)]
		frac.front[i,j] <- mean(ys.curr > 0 ,na.rm=T)
		median.dist.centr[i,j] <- median(sqrt(xs.curr^2 + ys.curr^2),na.rm=T)
		median.abs.x[i,j] <- median(abs(xs.curr),na.rm=T)
		median.abs.y[i,j] <- median(abs(ys.curr),na.rm=T)
	}
}

#consistency across days compared to randomized null
dat <- median.abs.x
#dat <- matrix(runif(dim(dat)[1]*dim(dat)[2]),nrow=dim(dat)[1],ncol=dim(dat)[2])
n.rands <- 100000
cons.dat <- mean(apply(dat,1,var,na.rm=T))
cons.rand <- rep(NA,n.rands)
for(i in 1:n.rands){
	dat.rand <- dat
	#shuffle non-na entries for each day
	for(j in 1:dim(dat)[2]){
		non.nas <- which(!is.na(dat[,j]))
		shuff <- sample(non.nas)
		dat.rand[non.nas,j] <- dat[shuff,j]
	}
	cons.rand[i] <- mean(apply(dat.rand,1,var,na.rm=T))
}

p <- mean(cons.rand < cons.dat)

#Individual consistency in fraction of time spent at front across days?
#Yes. cons.dat = .017, mean(cons.rand) = .021, p = .0002
#Using heading from 10m into the past: NO! cons.dat = .0108, mean(cons.rand)=.0109, p = .503

#Individual consistency in median distance from center?
#Yes. cons.dat = 4.15, mean(cons.rand) = 5.04, p = .00016
#Using heading from 10m into past: Yes. cons.dat = 3.93, mean(cons.rand) = 4.61, p = .00048

#Individual consistency in median abs distance to side?
#Yes. cons.dat = 2.06, mean(cons.rand) = 2.40, p = .0053
#Using heading from 10m into past: Yes. cons.dat=1.05, mean(cons.rand) = 1.33, p < .00001

#Individual consistency in median abs distance front/back?
#Yes. cons.dat = 1.78, mean(cons.rand) = 2.21, p < .001

#Individual consistency in fraction of time spent to the right of centroid across days?
#No evidence. cons.dat=.008, mean(cons.rand)==.0084, p = .15

#PLOT FRONTNESS FRACTIONS
par(bg='white',col='black',col.axis='black',col.lab='black',col.main='black')
dat <- median.abs.x
n.inds <- dim(dat)[1]
cols = rep('black',n.inds)
cols[which(meerkat.ids$dom & meerkat.ids$sex=='F')] <- '#FF0000'
cols[which(meerkat.ids$dom & meerkat.ids$sex=='M')] <- '#0000FF'
cols[which(!meerkat.ids$dom & meerkat.ids$sex=='F')] <- '#FFAA66'
cols[which(!meerkat.ids$dom & meerkat.ids$se=='M')] <- '#66AAFF'
meds <- apply(dat,1,median,na.rm=T)
ranks <- rev(order(meds))
upper.box <- apply(dat,1,function(x){return(quantile(x,0.75,na.rm=T))})
lower.box <- apply(dat,1,function(x){return(quantile(x,0.25,na.rm=T))})
upper.line <- apply(dat,1,function(x){return(quantile(x,0.975,na.rm=T))})
lower.line <- apply(dat,1,function(x){return(quantile(x,0.025,na.rm=T))})
plot(NULL,xlim=c(0,length(ranks)+1),ylim=c(1,7),xlab='',ylab='Median distance to side (m)',xaxt='n')
axis(1,at=1:n.inds,labels=meerkat.ids$id[ranks],las=2)
abline(h=50,lty=2)
arrows(1:n.inds,meds[ranks],1:n.inds,upper.line[ranks],length=0.1,lwd=2,angle=90,col='black')
arrows(1:n.inds,meds[ranks],1:n.inds,lower.line[ranks],length=0.1,lwd=2,angle=90,col='black')
for(i in 1:dim(dat)[1]){
	polygon(c(i-.25,i+.25,i+.25,i-.25),c(lower.box[ranks[i]],lower.box[ranks[i]],upper.box[ranks[i]],upper.box[ranks[i]]),col=cols[ranks[i]],lwd=2,border='black')
}
arrows((1:n.inds)-.25,meds[ranks],(1:n.inds)+.25,meds[ranks],length=0,lwd=3,col='black')


#PLOT FRAC FRONT VS DAY
plot(NULL,xlim=c(0,dim(frac.front)[2]+1),ylim=c(0,100),xlab='day',ylab='% time front of centroid')
abline(h=50,lty=1)
abline(v=15.5,lty=2)
abline(v=10.5,lty=2)
types <- c(1,1,1,1,2,3,6)
for(i in 1:n.inds){
	idxs <- which(!is.nan(frac.front[i,]))
	lines(idxs,frac.front[i,idxs]*100,col=cols[i],lwd=2,lty=types[i])
	points(frac.front[i,]*100,col=cols[i],pch=19,cex=1)	
}

