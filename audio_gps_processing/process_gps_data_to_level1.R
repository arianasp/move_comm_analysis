#interpole to fill small gaps in gps data to generate "level 1" dataset

#level0 dataset
level0.filename <- '~/Dropbox/meerkats/data/Kalahari2017/RDATA/gps_level0.RData'
max.gap.len <- 5 #maximum gap length to interpolate linearly (in sec)
max.speed <- 10 #maximum speed to consider reasonable (m/s)

#load in level 0 data
load(level0.filename)

#remove unrealistic locations (replace with NAs)
xs[which(xs<0)] <- NA
ys[which(xs<0)] <- NA
xs[which(ys<0)] <- NA
ys[which(ys<0)] <- NA

#remove cases where no x or no y coordinate exists
xs[which(is.na(ys))] <- NA
ys[which(is.na(xs))] <- NA

#remove unrealistic speeds (replace with NAs)
n.times <- dim(xs)[2]
n.inds <- dim(xs)[1]
speeds <- sqrt( ( xs[,2:n.times] - xs[,1:(n.times-1)] )^2 + ( ys[,2:n.times] - ys[,1:(n.times-1)] )^2 )
speeds <- cbind(rep(NA,n.inds),speeds)
xs[which(speeds>max.speed)] <- NA
ys[which(speeds>max.speed)] <- NA

#new x and y matrices
xs.new <- xs
ys.new <- ys

for(i in 1:dim(xs)[1]){
	print(i)
	for(d in 1:(length(day.start.idxs)-1)){
		xs.curr <- xs[i,day.start.idxs[d]:(day.start.idxs[d+1]-1)]
		ys.curr <- ys[i,day.start.idxs[d]:(day.start.idxs[d+1]-1)]
		tracked <- !is.na(xs.curr)

		idxs.tracked <- which(tracked)

		diffs <- diff(idxs.tracked)

		holes <- which(diffs > 1 & diffs <= (max.gap.len+1)) #find holes to fill

		if(length(holes > 0)){
			for(j in 1:length(holes)){
				t0 <- idxs.tracked[holes[j]] #index before gap
				tf <- idxs.tracked[holes[j]+(diffs[holes[j]]-1)] #index after gap
				x0 <- xs.curr[t0]
				y0 <- ys.curr[t0]
				xf <- xs.curr[tf]
				yf <- ys.curr[tf]
				x.vals <- seq(x0,xf,length.out=tf-t0+1)
				y.vals <- seq(y0,yf,length.out=tf-t0+1)
				xs.new[i,seq(day.start.idxs[d]+t0-1,day.start.idxs[d]+tf-1,1)] <- x.vals
				ys.new[i,seq(day.start.idxs[d]+t0-1,day.start.idxs[d]+tf-1,1)] <- y.vals
			}
		}
	}
}

xs <- xs.new
ys <- ys.new

save(list=c('xs','ys','times','day.start.idxs'),file='/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/xy_level1.RData')