#Heat map of locations used on n days (color represents n)

#LIBRARIES
library(raster)
library(dismo)

#SOURCE FUNCS
source('/Users/astrandb/Dropbox/hyenas/Hyena_pilot/Code/hyena_functions.R')

#LOAD DATA
load('/Users/astrandb/Dropbox/hyenas/HyenaGPS/RDATA/hyena_xy_level1.RData')
load('/Users/astrandb/Dropbox/hyenas/HyenaGPS/RDATA/hyena_timestamps.RData')
load('/Users/astrandb/Dropbox/hyenas/HyenaGPS/RDATA/hyena_ids.RData')

known.locs <- read.csv('/Users/astrandb/Dropbox/hyenas/HyenaGPS/CSV/hyena_isolate_dens.csv',stringsAsFactors=F)
eastsNorths <- latlon.to.utm(cbind(known.locs$lon,known.locs$lat),southern_hemisphere=T,utm.zone=36)
known.locs$east <- eastsNorths[,1]
known.locs$north <- eastsNorths[,2]

#parameters
bin.size <- 20 #size in meters of raster bins
#hours.to.include <- c(seq(0:6),seq(19:23)) #night
hours.to.include <- seq(7,18)

#get extent of x and y (easts and norths)
xmin <- quantile(xs,0.005,na.rm=T)
xmax <- quantile(xs,0.995,na.rm=T)
ymin <- quantile(ys,0.005,na.rm=T)
ymax <- quantile(ys,0.995,na.rm=T)

ext <- extent(xmin,xmax,ymin,ymax)

#get x and y bins for raster
x.bins <- seq(xmin,xmax,by=bin.size)
y.bins <- seq(ymin,ymax,by=bin.size)

#create a raster
rast <- raster()
extent(rast) <- ext
res(rast) <- bin.size
crs(rast) <- '+proj=utm +zone=36 +south +ellps=WGS84'

rast.dat <- data.frame()

#create a raster stack of each day and hyena
#also store information in rast.dat (hyena id and date of each raster in the raster stack)
rast.all <- stack()
dates <- date(timestamps+3*60*60) #correct time zone
dates.all <- unique(dates)
hours <- hour(timestamps+3*60*60) #hour of the day
for(i in 1:nrow(hyena.ids)){
	print(i)
	for(d in 1:length(dates.all)){
		idxs <- which(dates == dates.all[d] & hours %in% hours.to.include)
  		curr.x <- xs[i,idxs]
  		curr.y <- ys[i,idxs]
  
  		if(sum(!is.na(curr.x))>0){
			new.row <- data.frame(i=i,date=dates.all[d])
			rast.dat <- rbind(rast.dat,new.row)
			r0 <- rasterize(cbind(curr.x,curr.y),rast,fun=function(x,...){length(x)})
			rast.all <- stack(rast.all,r0)
  		}
	}
}

#create a heat map for each hyena - how many days was each raster bin used by that ind?

ind.rasts <- stack()
for(ind in 1:nrow(hyena.ids)){
	new.rast <- raster(res=res(rast),ext=extent(rast))

	for(i in which(rast.dat$i==ind)){
		new.rast <- mosaic(new.rast,rast.all[[i]]>0,fun='sum')
	}
	ind.rasts <- stack(ind.rasts,new.rast)
}

locs.to.plot <- known.locs[which(known.locs$name %in% c('DAVE D','RES M D1','RBEND D')),]
ds <- grep(' D', known.locs$name)
nds <- grep(' ND',known.locs$name)
cds <- grep(' CD',known.locs$name)
cols.dens <- c('red','yellow','green')
known.locs$col <- 'black'
known.locs$col[ds] <- 'red'
known.locs$col[nds] <- 'green'
known.locs$col[cds] <- 'yellow'

quartz(width=10,height=8)
par(mfrow=c(2,3),mar=c(0,0,0,0))
for(i in 1:nrow(hyena.ids)){
	image(ind.rasts[[i]],col=viridis(256),asp=1,xaxt='n',yaxt='n',zlim=c(0,29))
	#points(known.locs$east,known.locs$north,col=known.locs$col,cex=0.5)
	#points(locs.to.plot$east,locs.to.plot$north,col='yellow',cex=0.5)
	text(locs.to.plot$east,locs.to.plot$north,col='yellow',locs.to.plot$name,cex=0.5,pos=4)
	textx <- extent(rast)[1]+(extent(rast)[2]-extent(rast)[1])/10
	texty <- extent(rast)[3]+(extent(rast)[4]-extent(rast)[3])/10
	text(textx,texty,hyena.ids$name[i],col='white')	
}

#make a meta map (stack the stacks, take only locations used at least 3 times)
meta.rast <- raster(res=res(rast),ext=extent(rast))
for(i in 1:nrow(hyena.ids)){
	curr.rast <- ind.rasts[[i]] > 3
	meta.rast <- mosaic(meta.rast,curr.rast,fun='sum')
}
ind.space.rast <- meta.rast==1
shared.space.rast <- meta.rast > 1
image(shared.space.rast*meta.rast,col=viridis(256),asp=1,xaxt='n',yaxt='n',zlim=c(0,5))
#points(known.locs$east,known.locs$north,col=known.locs$col,cex=0.5)
text(locs.to.plot$east,locs.to.plot$north,col='yellow',locs.to.plot$name,cex=0.5,pos=4)
textx <- extent(rast)[1]+(extent(rast)[2]-extent(rast)[1])/10
texty <- extent(rast)[3]+(extent(rast)[4]-extent(rast)[3])/10
text(textx,texty,hyena.ids$name[i],col='white')	

#plot on google earth imagery
r2 <- projectRaster(rast,crs='+proj=utm +zone=36 +south +ellps=WGS84')
g <- gmap(r2,type='satellite',lonlat=F,rgb=T)
g2 <- projectRaster(g,crs='+proj=utm +zone=36 +south +ellps=WGS84',res=c(20,20))
g3 <- crop(g2,rast) 
extent(g3) <- extent(rast)

#Make plots on maps!
quartz(width=10,height=8)
par(mfrow=c(2,3),mar=c(0,0,0,0))
for(i in 1:nrow(hyena.ids)){
	maxi <- maxValue(ind.rasts[[i]])
	mini <- minValue(ind.rasts[[i]])
	curr.rast <- mosaic(ind.rasts[[i]]*500/(maxi-mini),g3,fun='sum')
	curr.rast <- mosaic(ind.rasts[[i]]*20,g3,fun='sum')
	plotRGB(curr.rast,scale=680)
	#points(known.locs$east,known.locs$north,col=known.locs$col,cex=0.5)
	#points(locs.to.plot$east,locs.to.plot$north,col='yellow',cex=0.5)
	text(locs.to.plot$east,locs.to.plot$north,col='yellow',locs.to.plot$name,cex=0.5,pos=4)
	textx <- extent(rast)[1]+(extent(rast)[2]-extent(rast)[1])/10
	texty <- extent(rast)[3]+(extent(rast)[4]-extent(rast)[3])/10
	text(textx,texty,hyena.ids$name[i],col='white')	
}
maxi <- maxValue(meta.rast)
mini <- minValue(meta.rast)
curr.rast <- mosaic(meta.rast*300/(maxi-mini),g3,fun='sum')
plotRGB(curr.rast,scale=680)
#points(known.locs$east,known.locs$north,col=known.locs$col,cex=0.5)
text(locs.to.plot$east,locs.to.plot$north,col='yellow',locs.to.plot$name,cex=0.5,pos=4)
text(textx,texty,'SHARED',col='white')	
