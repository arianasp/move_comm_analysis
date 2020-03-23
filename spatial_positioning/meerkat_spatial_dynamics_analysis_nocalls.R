#meerkat spatial dynamics analysis

#do meerkats move toward or away from each other as a function of distance?

#make a meerkat 'force map' aligned with the group heading

#bring in calls on days when we have them - does the curve shift with call rate?

#make a map of call rate vs. normalize spatial position in group

library(pracma)
library(fields)

#SOURCE
source('/Users/astrandb/Dropbox/meerkats/code/general/meerkat_functions.R')
source('/Users/astrandb/Dropbox/code_ari/baboons_code/step_selection/step_selection_funcs.R')

#LOAD DATA
load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/positions_rel_to_centroid.RData')
load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/xy_level1.RData')
load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/ids.RData')
load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/babysit_periods.RData')
load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/calls_all.RData')

#REMOVE THE BABYSITTERS
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

#EXTRACT DATA
#for each pair of individuals and at each moment in time (moving forward in step sec intervals)
#calculate how much time passes before the individual moves R m away
step <- 10
R <- 5
xbins <- ybins <- c(-50,-20,-10,-5,-2,-1,0,1,2,5,10,20,50)
max.bin <- max(xbins)
leave.times <- leave.direcs <- leave.twd.away.vels <- dists <- dxs.rel <- dys.rel <- dxs.abs <- dys.abs <- is <- js <- dxs.past <- dys.past <- angs.past <- array(NA,dim=c(length(xbins),length(ybins),n.times/5))
counts <- array(0,dim=c(length(xbins),length(ybins)))

for(i in 1:n.inds){ #i is the mover
	print(i)
	for(j in (1:n.inds)[-i]){ #j is the reference individual
		print(j)
		for(d in 1:(length(day.start.idxs)-1)){
			print(d)
			trange <- seq(day.start.idxs[d],day.start.idxs[d+1]-1) #times associated with each day
			xs.curr <- xs[,trange] #x values for that day
			ys.curr <- ys[,trange] #y avlues for that day
			xs.rel.curr <- xs.rel[,trange] #relative x values for that day
			ys.rel.curr <- ys.rel[,trange] #relative y values for that day
			for(t in seq(1,length(trange),step)){
			
				#get initial x and y of each individual relative to the group reference frame
				xi.rel <- xs.rel.curr[i,t]
				yi.rel <- ys.rel.curr[i,t]
				xj.rel <- xs.rel.curr[j,t]
				yj.rel <- ys.rel.curr[j,t]
				
				#initial vector pointing from focal to reference individual - group reference frame
				dx.rel <- xj.rel - xi.rel 
				dy.rel <- yj.rel - yi.rel
				
				#initial x and y of focal individual in static reference frame
				xi <- xs.curr[i,t]
				yi <- ys.curr[i,t]
				xj <- xs.curr[j,t]
				yj <- ys.curr[j,t]
				
				#initial vector point from focal to reference individual - static reference frame
				dx.abs <- xj - xi
				dy.abs <- yj - yi
				
				#initialize time of leaving from a radius R to NA
				time <- NA
				
				if(!is.na(dx.rel)){
					if(abs(dx.rel) < max.bin & abs(dy.rel) < max.bin){
						xbin.idx <- which(xbins > dx.rel)[1]
						ybin.idx <- which(ybins > dy.rel)[1]
						
						#find the time when the focal leaves a radius R from its current location
						found <- F
						fwd.idx <- 1
						while(!found){
						
							#quit if you reach the end of the day, otherwise continue
							if(((t+fwd.idx) > dim(xs.curr)[2])){
								found <- T
								time <- NA
							} else{
							
								#get distance of current focal location (at t+fwd.idx) from initial location (at t)
								dist.curr <- sqrt( (xi - xs.curr[i,t+fwd.idx])^2 + (yi - ys.curr[i,t+fwd.idx])^2 )
								
								#if you run into an NA, quit
								if(is.na(dist.curr)){
									found <- T
									data.found <- F
									time <- NA
								} else{
									
									#if the distance exceeds R, use this time and location for the rest of the calculations
									if(dist.curr > R){
									
										#time it took to leave the circle
										time <- fwd.idx 
										
										#leaving velocity vector of focal (in m/s) - static reference frame
										leave.dx <- (xs.curr[i,t+fwd.idx] - xs.curr[i,t]) / time
										leave.dy <- (ys.curr[i,t+fwd.idx] - ys.curr[i,t]) / time
										
										#vector pointing from focal (mover) initial location to reference individiual initial location - static reference frame
										dx.foc.ref <- xs.curr[j,t] - xs.curr[i,t]
										dy.foc.ref <- ys.curr[j,t] - ys.curr[i,t]
										
										#vector pointing from focal (mover) to reference individaul (normalized) - static reference frame
										dx.norm <- dx.foc.ref / sqrt(dx.foc.ref^2+dy.foc.ref^2)
										dy.norm <- dy.foc.ref / sqrt(dx.foc.ref^2+dy.foc.ref^2)
										
										#normalized leaving vector of focal - static reference frame
										leave.dx.norm <- leave.dx / sqrt(leave.dx^2 + leave.dy^2)
										leave.dy.norm <- leave.dy / sqrt(leave.dx^2 + leave.dy^2)
										
										#directionality of leaving (from - 1 = away from ref ind to 1 = toward ref ind) - static reference frame
										direc <- leave.dx.norm*dx.norm + leave.dy.norm*dy.norm
										
										#directionality * velocity of leaving - static reference frame
										twd.away <- leave.dx*dx.norm + leave.dy*dy.norm
										
										#set found flag true
										found <- T
										data.found <- T #flag that indicates a non-NA value will be returned
										
										#store data
										counts[xbin.idx,ybin.idx] <- counts[xbin.idx,ybin.idx]+1 #events in each bin
										next.spot <- counts[xbin.idx,ybin.idx] #where to store in the arrays
										leave.times[xbin.idx,ybin.idx,next.spot] <- time #store leaving time
										leave.direcs[xbin.idx,ybin.idx,next.spot] <- direc #store leaving direction
										leave.twd.away.vels[xbin.idx,ybin.idx,next.spot] <- twd.away #store leaving toward/away velocity
										dists[xbin.idx,ybin.idx,next.spot] <- sqrt(dx.foc.ref^2+dy.foc.ref^2) #store initial distance between focal and reference individual
										dxs.rel[xbin.idx,ybin.idx,next.spot] <- dx.rel
										dys.rel[xbin.idx,ybin.idx,next.spot] <- dy.rel
										dxs.abs[xbin.idx,ybin.idx,next.spot] <- dx.abs
										dys.abs[xbin.idx,ybin.idx,next.spot] <- dy.abs
										is[xbin.idx,ybin.idx,next.spot] <- i
										js[xbin.idx,ybin.idx,next.spot] <- j
										
									} else{
	
										#increment time
										fwd.idx <- fwd.idx+1
										
									}
								}
							}
						}
						
						#get previous location (R m into the past) to get initial heading
						bwd.idx <- 1
						found <- F
						if(data.found){
								while(!found){
						
								#quit if you reach the beginning of the day, otherwise continue
								if(((t-bwd.idx) < 1)){
									found <- T
									time <- NA
								} else{
							
									#get distance of current focal location (at t+fwd.idx) from initial location (at t)
									dist.curr <- sqrt( (xi - xs.curr[i,t-bwd.idx])^2 + (yi - ys.curr[i,t-bwd.idx])^2 )
								
									#if you run into an NA, quit
									if(is.na(dist.curr)){
										found <- T
										time <- NA
									} else{
									
										#if the distance exceeds R, use this time and location for the rest of the calculations
										if(dist.curr > R){
									
											#time it took to leave the circle
											time <- bwd.idx 
										
											#initial velocity vector of focal (in m/s)
											init.dx <- -(xs.curr[i,t-bwd.idx] - xs.curr[i,t]) / time
											init.dy <- -(ys.curr[i,t-bwd.idx] - ys.curr[i,t]) / time
										
											#normalized initial heading vector of focal
											init.dx.norm <- init.dx / sqrt(init.dx^2 + init.dy^2)
											init.dy.norm <- init.dy / sqrt(init.dx^2 + init.dy^2)
										
											#set found flag true
											found <- T
										
											#store data
											dxs.past[xbin.idx,ybin.idx,next.spot] <- init.dx
											dys.past[xbin.idx,ybin.idx,next.spot] <- init.dy
											angs.past[xbin.idx,ybin.idx,next.spot] <- atan2(init.dy,init.dx)
											
										} else{
	
											#increment time
											bwd.idx <- bwd.idx+1
										
										}
									}
								}
							}
						}
					}
				}
			}	
		}
	}
}

save(file='/Users/astrandb/Desktop/leave_times3.RData',list=ls())


#MAKE PLOTS


#get data into lists to match up with distances 
all.dists <- c(dists)
all.leave.times <- c(leave.times)
all.leave.direcs <- c(leave.direcs)
all.twd.away.vel <- c(leave.twd.away.vels)*60 #convert to m/min
all.speeds <- 5 / all.leave.times * 60 #speed of leaving radius = 5
all.is <- c(is)
all.js <- c(js)
all.dxs.past <- c(dxs.past)
all.dys.past <- c(dys.past)
all.angs.past <- c(angs.past)
all.dxs.rel <- c(dxs.rel)
all.dys.rel <- c(dys.rel)
all.dxs.abs <- c(dxs.abs)
all.dys.abs <- c(dys.abs)
all.angs.rel <- atan2(all.dys.rel,all.dxs.rel)
all.angs.abs <- atan2(all.dys.abs,all.dxs.abs)
all.angs.prev <- atan2(all.dys.past,all.dxs.past)
all.angle.changes <- all.angs.abs - all.angs.prev
all.angle.changes[which(all.angle.changes > pi)] <- all.angle.changes[which(all.angle.changes > pi)] - 2*pi
all.angle.changes[which(all.angle.changes < -pi)] <- all.angle.changes[which(all.angle.changes < -pi)] + 2*pi


y <- all.leave.direcs #stat to use


#MAKE 2D PLOTS OF DISTANCE VS LEAVE TIME / DIRECTION / ETC
dist.bins <- quantile(all.dists,seq(0,.95,.05),na.rm=T) #distance bins
ylab = '% approach'

medians <- means <- fracs <-  q25 <- q05 <- q75 <- q95 <- counts <- rep(NA,length(dist.bins)-1)
for(i in 1:(length(dist.bins)-1)){
	idxs <- which(all.dists >= dist.bins[i] & all.dists < dist.bins[i+1])
	medians[i] <- median(y[idxs],na.rm=T)
	means[i] <- mean(y[idxs],na.rm=T)
	q05[i] <- quantile(y[idxs],.05,na.rm=T)
	q25[i] <- quantile(y[idxs],.25,na.rm=T)
	q75[i] <- quantile(y[idxs],.75,na.rm=T)
	q95[i] <- quantile(y[idxs],.95,na.rm=T)
	counts[i] <- length(idxs)
	fracs[i] <- mean(y[idxs]>0,na.rm=T)
}
par(bg='black',col='white',col.axis='white',col.lab='white',col.main='white')

#make the plot
x <- (dist.bins[2:length(dist.bins)] + dist.bins[1:(length(dist.bins)-1)])/2
y.plot <- fracs*100
plot(NULL,xlab='initial distance between meerakts (m)',ylab=ylab,xlim=c(0,max(dist.bins)),ylim=c(50,65))
#polygon(c(x,rev(x)),c(q05,rev(q95)),col='#00000022',border=NA)
#polygon(c(x,rev(x)),c(q25,rev(q75)),col='#00000022',border=NA)
#lines(x,means,lwd=3)
arrows(dist.bins[1:(length(dist.bins)-1)],y.plot,dist.bins[2:length(dist.bins)],y.plot,length=0,lwd=3,col='white')
abline(h=50,lty=2)



#MAKE 3D PLOTS (color)

y <- all.leave.direcs #stat to use

mean.leave.times <- apply(leave.times,c(1,2),mean,na.rm=T)
median.leave.times <- apply(leave.times,c(1,2),median,na.rm=T)
mean.twd.away <- apply(leave.direcs,c(1,2),mean,na.rm=T)
mean.twd.away.vel <- apply(leave.twd.away.vels,c(1,2),mean,na.rm=T)

cols.palette <- colorRampPalette(c('red','white','blue'))
cols <- cols.palette(256)


z <- mean.twd.away
image.plot(z[2:(length(xbins)-1),2:(length(xbins)-1)],zlim=c(-1,1),col=cols)
abline(h=(length(xbins)-1)/2,col='white')
abline(v=(length(xbins)-1)/2,col='white')

#circular heat maps


ang.bins <- seq(-pi,pi,pi/8)
dist.bins <- quantile(all.dists,seq(0,.95,.05),na.rm=T)
medians.c <- means.c <- counts.c <- fracs.c <- array(NA,dim=c(length(dist.bins)-1,length(ang.bins)-1))
approachee <- 7
y <- all.leave.direcs
for(i in 1:(length(dist.bins)-1)){
	for(j in 1:(length(ang.bins)-1)){
		#idxs <- which((all.dists >= dist.bins[i]) & (all.dists < dist.bins[i+1]) & (all.angs.rel >= ang.bins[j]) & (all.angs.rel < ang.bins[j+1]) & (all.js==approachee)) #using angle relative to group heading
		idxs <- which((all.dists >= dist.bins[i]) & (all.dists < dist.bins[i+1]) & (all.angle.changes >= ang.bins[j]) & (all.angle.changes < ang.bins[j+1]) & (all.js==approachee)) #using angle relative to previous individual heading
		medians.c[i,j] <- median(y[idxs],na.rm=T)
		means.c[i,j] <- mean(y[idxs],na.rm=T)
		counts.c[i,j] <- length(idxs)
		fracs.c[i,j] <- mean(y[idxs]>0)
	}
}

#make a circular plot
z <- fracs.c
zlim <- c(0,1)
col.palette <- colorRampPalette(c('red','white','blue'))
cols <- col.palette(1024)
plot(NULL,xlim=c(-max(dist.bins),max(dist.bins)),ylim=c(-max(dist.bins),max(dist.bins)),xlab='x relative (m)',ylab='y relative (m)',asp=1)
for(i in 1:(length(dist.bins)-1)){
	for(j in 1:(length(ang.bins)-1)){
		dat <- z[i,j]
		col <- cols[round((dat - zlim[1])/(zlim[2] - zlim[1])*length(cols))]
		r0 <- dist.bins[i]
		rf <- dist.bins[i+1]
		ang0 <- ang.bins[j]
		angf <- ang.bins[j+1]
		d.angs <- seq(ang0,angf,(angf - ang0)/100)
		x1 <- r0*cos(ang0)
		x2 <- rf*cos(d.angs)
		x3 <- r0*cos(angf)
		x4 <- r0*cos(rev(d.angs))
		y1 <- r0*sin(ang0)
		y2 <- rf*sin(d.angs)
		y3 <- r0*sin(angf)
		y4 <- r0*sin(rev(d.angs))
		polygon(c(x1,x2,x3,x4),c(y1,y2,y3,y4),col=col,border=NA)
	}
}

plot(NULL,xlim=c(min(dist.bins),max(dist.bins)),ylim=c(min(z,na.rm=T),max(z,na.rm=T)))
cols=rainbow(length(ang.bins))
for(i in 1:(length(ang.bins)-1)){
	lines(dist.bins[2:length(dist.bins)],z[,i],col=cols[i])
}

#make plots by individual
y <- all.twd.away.vel
dist.bins <- quantile(all.dists,seq(0,1,.2),na.rm=T)
max.abs.ang <- pi/8 #for getting only the sides
means <- medians <- fracs <- counts <- array(NA,dim=c(n.inds,length(dist.bins)-1))
for(j in 1:n.inds){
	for(k in 1:(length(dist.bins)-1)){
		idxs <- which((all.js==j) & (all.dists > dist.bins[k]) & (all.dists <= dist.bins[k+1]) & (abs(cos(all.angs.rel)) > cos(max.abs.ang)))
		dat <- y[idxs]
		means[j,k] <- mean(dat,na.rm=T)
		medians[j,k] <- median(dat,na.rm=T)
		fracs[j,k] <- mean(dat>0,na.rm=T)
		counts[j,k] <- length(idxs)
	}
}
par(bg='black',col='white',col.axis='white',col.lab='white',col.main='white')
y.plot <- fracs*100
x.plot <- (dist.bins[2:length(dist.bins)] + dist.bins[1:(length(dist.bins)-1)])/2
cols = c('red','blue','#FFAA66','#66AAFF','#66AAFF','#66AAFF','#66AAFF')
plot(NULL,xlim=c(0,max(x.plot)+5),ylim=c(min(y.plot),max(y.plot)),ylab='% approach',xlab='initial distance between meerkats (m)',bg='black')
for(i in 1:n.inds){
	lines(x.plot,y.plot[i,],col=cols[i],lwd=3)
	text(max(x.plot)+3.2,y.plot[i,dim(y.plot)[2]],meerkat.ids$id[i],col=cols[i])
}

x <- (dist.bins[2:length(dist.bins)] + dist.bins[1:(length(dist.bins)-1)])/2
plot(NULL,xlim=c(0,max(x)),ylim=c(min(fracs)*100,max(fracs)*100),xlab='initial distance between meerkats (m)',ylab='% approach')
cols = c('red','blue','#FFAA66','#66AAFF','#66AAFF','#66AAFF','#66AAFF')
for(j in 1:n.inds){
	lines(x,fracs[j,]*100,col=cols[j],lwd=2)
	#points(x,fracs[j,]*100,col=cols[j],pch=19)
}
abline(h=50,lty=2)