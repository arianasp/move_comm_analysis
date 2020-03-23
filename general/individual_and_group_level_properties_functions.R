#This file contains a set of functions useful for analyzing collective movement data. 
#The functions compute both individual and group-level properties of animal groups.
#I have tried to make the functions human-readable so that you may understand how they work, 
#and have tried to keep the implementation as simple as possible for this purpose as well.

#NOTE: Note: use the 'fields' library to produce prettier heat map plots! Uncomment-out the line below 
#(first you must download it using install.packages('fields')):
#library(fields)

#-----------------INDIVIDUAL-LEVEL PROPERTIES-----------------------

#FUNCTION: get.speed.over.time
#This function gets the instantaneous speed of an individual over time. 
#The speed is computed as s[t] = sqrt ( (x.i[t+1] - x.i[t] )^2 + ( y.i[t+1] - y.i[t] )^2 ) / dt
#Inputs:
#	x.i: x coordinates of the individual (a vector whose length is the number of timesteps) or of the group centroid
#	y.i: y coordinates of the individual (a vector whose length is the number of timesteps) or of the group centroid
#	frame.rate: the number of frames per second in the input data (defaults to 20)
#Outputs:
#	speeds: time series of the speed of individual
get.speed.over.time <- function(x.i,y.i,frame.rate=20){
	
	#get the change in x and y coordinates for each time step
	dx <- diff(x.i) 
	dy <- diff(y.i)
	
	#append an NA to the end of the distance vectors to keep them the same length as the input data
	dx <- c(dx,NA)
	dy <- c(dy,NA)
	
	#get change in time
	dt <- 1 / frame.rate
	
	#get instantaneous speed
	speeds <- sqrt( dx^2 + dy^2 ) / dt
	
	#output speeds
	return(speeds)
	
}

#FUNCTION: get.nearest.neighbor.distance
#Compute the distance to an individual's nearest neighbor over time
#This function gets the distance to an individual i's nearest neighbor over time
#Inputs:
#	xs: N x n.times matrix giving x positions of all N individuals over n.times timesteps
#	ys: N x n.times matrix giving y positions of all N individuals over n.times timesteps
#	i: the focal individual (whose nearest neighbor distance will be computed)
#Outputs:
#	nn.dists: vector of length n.times giving the distance to the focal individual's nearest neighbor at each time point
get.nearest.neighbor.distance <- function(xs,ys,i){
	
	#get the position of individual i over time
	xs.i <- xs[i,]
	ys.i <- ys[i,]
	
	#number of timesteps and number of individuals
	N <- nrow(xs)
	n.times <- ncol(xs)
	
	#convert these positions to matrices (which we will subtract in the next step to get distances to all individuals)
	xs.i.mat <- matrix(rep(xs.i,each=N),nrow=N,ncol=n.times)
	ys.i.mat <- matrix(rep(ys.i,each=N),nrow=N,ncol=n.times)
	
	#subtract to get distances in x and y dimensions between focal individual and all other individuals
	dx.mat <- xs - xs.i.mat
	dy.mat <- ys - ys.i.mat
	
	#get total distance between focal and all other individuals
	dists.mat <- sqrt( dx.mat^2 + dy.mat^2 )
	
	#distance between a focal individual and itself will be 0, so set this to NA
	dists.mat[i,] <- NA
	
	#find the minimum distance at each time (column) - this is the nearest neighbor distance
	nn.dists <- apply(dists.mat,MARGIN=2,FUN=min,na.rm=T)
	
}

#FUNCTION: get.local.density
#This function gives the 'local density' or the number of individuals within a radius R of a focal individual
#Inputs:
#	xs: N x n.times matrix giving x positions of all N individuals over n.times timesteps
#	ys: N x n.times matrix giving y positions of all N individuals over n.times timesteps
#	i: the focal individual (for whom the local density will be computed)
#	R: the radius over which to compute the local density (in same units as xs and ys)
#Outputs:
#	local.density: vector of length n.times giving local density around focal individual (number of individuals within a radius R) at each time
get.local.density <- function(xs,ys,i,R){
	
	#A quick note: this function starts out exactly the same as the one above 
	#(i.e. getting a matrix of distances from focal individual to all others in the group)
	#This perhaps means that it would be a good idea to create a separate function to compute this matrix,
	#and then use it within both of these functions. I haven't done this here for the sake of clarity of the code,
	#but be aware that functions calling other functions is a regular and useful thing!
	
	#get the position of individual i over time
	xs.i <- xs[i,]
	ys.i <- ys[i,]
	
	#number of timesteps and number of individuals
	N <- nrow(xs)
	n.times <- ncol(xs)
	
	#convert these positions to matrices (which we will subtract in the next step to get distances to all individuals)
	xs.i.mat <- matrix(rep(xs.i,each=N),nrow=N,ncol=n.times)
	ys.i.mat <- matrix(rep(ys.i,each=N),nrow=N,ncol=n.times)
	
	#subtract to get distances in x and y dimensions between focal individual and all other individuals
	dx.mat <- xs - xs.i.mat
	dy.mat <- ys - ys.i.mat
	
	#get total distance between focal and all other individuals
	dists.mat <- sqrt( dx.mat^2 + dy.mat^2 )
	
	#distance between a focal individual and itself will be 0, so set this to NA
	dists.mat[i,] <- NA
	
	#find the number of individuals withiin a distance R at each time point (column)
	local.density <- colSums(dists.mat <= R,na.rm=T)
	
}

#FUNCTION: get.directedness
#Gets the directedness of a trajectory, a number which ranges from 0 to 1 where 1 is a straight path and 0 is a highly tortuous path
#The directedness is defined as the net displacement (distance along a straight-line path from point A to point B) divided by the
# path length (total distance traveled along the actual path from point A to point B)
#Inputs:
#	x.i: vector of x positions for an individual (or the group centroid) of length n.times
#	y.i: vector of y positions for an individual (or the group centroid) of length n.times
#	t.window: window of time to use for computing directedness (must be an even number) - the directedness for time t will be computed using position data from t - t.window/2 to t + t.window / 2
#NOTE: Keep in mind that the time window should be given as the number of time points of data to use, so the unit is frames not seconds!
#Outputs:
#	directedness: vector of directedness of the trajectory as a function of time
get.directedness <- function(x.i, y.i,t.window){
	
	#number of time points in the position vectors
	n.times <- length(x.i)
	
	#create a vector to hold directedness values
	directedness <- rep(NA,n.times)
	
	#this is called a "for loop" - everything inside gets executed for a range of values of the "index" variable, in this case "t"
	#in this case, i indexes over the time, so we are computing the directedness for each value of the time
	#we often start for loops at 1, but here we will actually start at t.window/2 + 1 because we cannot compute the directedness
	#for the first t.window / 2 time points, since this would require having data from before the beginning of the file.
	
	for(t in seq(t.window/2+1,n.times - t.window/2,1)){ #seq makes a sequence the runs from t.window/2+1 to n.times - t.window/2 by steps of 1
		
		#everything inside the curly brackets (these things {}) is part of the loop
		
		#get values of x and y positions for the individual within the current time window
		x.curr <- x.i[seq(t-t.window/2+1,t+t.window/2,1)]
		y.curr <- y.i[seq(t-t.window/2+1,t+t.window/2,1)]
		
		#get change in x and y positions during each time step
		dx.curr <- diff(x.curr)
		dy.curr <- diff(y.curr)
		
		#get distance traveled in each time step
		dists <- sqrt(dx.curr^2 + dy.curr^2)
		
		#get the total distance traveled ("path length") across the time window
		path.len <- sum(dists)
		
		#get the net displacement (distance from the position at beginning of the time window to position at end of time window)
		dx.tot <- x.curr[length(x.curr)] - x.curr[1] #displacement in x dimension
		dy.tot <- y.curr[length(y.curr)] - y.curr[1] #displacement in y dimension
		net.disp <- sqrt( dx.tot^2 + dy.tot^2 )
		
		#get the directedness at that time point, which is net.disp / path.len
		directedness.curr <- net.disp / path.len
		
		#store this directedness at position i in the directedness vector (since this is the directedness for the current time)
		directedness[t] <- directedness.curr
		
		#now the loop repeats with the next value of t!
		
	} 
	
	#at the end of the loop, output the vector 'directedness'
	return(directedness)
	
}


#---------------------GROUP LEVEL PROPERTIES---------------------------

#FUNCTION: get.group polarization
#Computes the polarization of the group at each time t. 
#The polarization is a measure of how aligned the group is, ranging from 0 (completely unaligned) to 1 (completely aligned)
#The polarization is defined by adding up (vector addition) all of the heading vectors of all individuals at a given moment in time, 
#taking the length of the resultant vector, and dividing this by the number of individuals, N. 
#Inputs:
#	xs: N x n.times matrix giving x coordinates of each individual over time
#	ys: N x n.times matrix giving y coordinates of each individual over time
#Output:
#	polarization: vector of length n.times giving the polarization over time 
#	(for the last point the heading is undefined, so the last element will be NA)

get.group.polarization <- function(xs,ys){
	
	#get displacements in x and y directions
	dx <- t(apply(xs,MARGIN=1,FUN=diff))
	dy <- t(apply(ys,MARGIN=1,FUN=diff))
	
	#get total displacement
	total.disp <- sqrt(dx^2 + dy^2)
	
	#get heading vectors (x and y components) by dividing displacement in each dimension by total displacement
	heads.x <- dx / total.disp
	heads.y <- dy / total.disp
	
	#add up heading vectors for each moment in time (column)
	heads.x.tot <- colSums(heads.x)
	heads.y.tot <- colSums(heads.y)
	
	#get length of the resultant vector
	resultant.vec.len <- sqrt(heads.x.tot^2 + heads.y.tot^2)
	
	#divide by the number of individuals (N) to get the polarization at each moment in time
	N <- nrow(xs)
	polarization <- resultant.vec.len / N
	
	#append an NA because the polarization is undefined for the last time index
	polarization <- c(polarization,NA)
	
	#output polarization over time
	return(polarization)
	
}

#FUNCTION: get.group.rotation
#This function computes the rotation of the group over time
#The rotation is a measure of how much the group rotates around its centroid, ranging from 0 (no rotation) to 1 (everyone rotating in the same direction)
#Inputs:
#	xs: N x n.times matrix giving x coordinates of each individual over time
#	ys: N x n.times matrix giving y coordinates of each individual over time
#Output:
#	rotation: vector of group rotation over time

get.group.rotation <- function(xs,ys){
	
	#number of individuals and times
	N <- nrow(xs)
	n.times <- ncol(xs)
	
	#first we need to get the headings (as in the last function):
	
	#get displacements in x and y directions
	dx <- t(apply(xs,MARGIN=1,FUN=diff))
	dy <- t(apply(ys,MARGIN=1,FUN=diff))
	
	#get total displacement
	total.disp <- sqrt(dx^2 + dy^2)
	
	#get heading vectors (x and y components) by dividing displacement in each dimension by total displacement
	heads.x <- dx / total.disp
	heads.y <- dy / total.disp
	
	#add a column of NAs for the last time point (when headings are undefined)
	heads.x <- cbind(heads.x,rep(NA,N))
	heads.y <- cbind(heads.y,rep(NA,N))
	
	#now we will need to compare these headings with the vector pointing from the group centroid to each individual
	
	#to do this, let's first compute the group centroid:
	centr.xs <- colMeans(xs)
	centr.ys <- colMeans(ys)
	
	#we will then make these into matrices to enable direct subtraction in the next step
	centr.x.mat <- matrix(rep(centr.xs,each=N),nrow=N,ncol=n.times)
	centr.y.mat <- matrix(rep(centr.ys,each=N),nrow=N,ncol=n.times)
	
	#now compute the direction vector from the centroid at time t to each individual's position at time t
	dxs.ind.centr <- xs - centr.x.mat
	dys.ind.centr <- ys - centr.y.mat
	
	#normalize these vectors so they have length 1 (divide x and y components by the distance from each individual to the centroid)
	dist.ind.centr <- sqrt( dxs.ind.centr^2 + dys.ind.centr^2)
	dxs.ind.centr.norm <- dxs.ind.centr / dist.ind.centr
	dys.ind.centr.norm <- dys.ind.centr / dist.ind.centr
	
	#take the cross product of these two vectors (and individual's heading and the direction from centroid to its position)
	cross.products <- dxs.ind.centr.norm*heads.y - dys.ind.centr.norm*heads.x
	
	#now compute the rotation as the sum of cross products across all individuals, divided by the number of individuals
	rotation <- colSums(cross.products) / N
	
}

#FUNCTION: get.mean.dyadic.distance
#This function computes the mean dyadic distance of the group over time.
#The mean dyadic distance is computed by taking the distance between every pair of individuals and then averaging these.
#Inputs:
#	xs: N x n.times matrix giving x coordinates of each individual over time
#	ys: N x n.times matrix giving y coordinates of each individual over time
#Outputs:
#	mean.dyadic.dist: vector of length n.times giving the mean dyadic distance as a function of time

get.mean.dyadic.dist <- function(xs,ys){
	
	#get the number of time points
	n.times <- ncol(xs)
	
	#create a vector to hold mean dyadic distance over time
	mean.dyad.dists <- rep(NA,n.times)
	
	#use a "for loop" to compute the mean dyadic distance at each point in time
	
	for(t in 1:n.times){ #for each time point...
		
		#get the x and y coordinates of all individuals at that time
		xs.t <- xs[,t]
		ys.t <- ys[,t]
		
		#get pairwise distances between every pair of individuals
		dyad.dists.t <- dist(cbind(xs.t,ys.t))
		
		#get the mean dyadic distance
		mean.dyad.dist.t <- mean(dyad.dists.t)
		
		#store this value in the output matrix
		mean.dyad.dists[t] <- mean.dyad.dist.t
		
	}
	
	#output the mean dyadic distance over time
	return(mean.dyad.dists)
	
}

#------------------PLOTTING FUNCTIONS----------------------
#FUNCTION: hist.2d
#This function makes a two-dimensional histogram plot with color representing probability
#Inputs:
#	x: vector of values corresponding to the first measure in the histogram (e.g. polarization)
#	y: vector of values correpsonding to the second measure in the histogram (e.g. rotation)
#	xbins: vector of bins for first measure (x) - defaults to seq(0,1,0.05)
#	ybins: vector of bins for second measure (y) - defaults to seq(0,1,0.05)
#	output_freqs: T or F, whether to output overall frequencies (if T) or probabilities of occurrence (if F) in the histogram - defaults to F
#	output_plot: T or F, whether to make a plot - defaults to T
#	xlab: label (must be a character string enclosed with quotes) for the x axis - defaults to no label
#	ylab: label (must be a character string enclosed with quotes) for the y axis - defaults to no label
#Outputs:
#	histo: matrix representing two-dimensional histogram of the data

hist.2d <- function(x,y,z = NULL, xbins=seq(0,1,0.05),ybins=seq(0,1,0.05),output_freqs=F,output_plot=T,xlab='',ylab='',xaxt=NULL,yaxt=NULL,axes=T,imageplot=T,zlim=NULL){
	nx <- length(xbins)
	ny <- length(ybins)
	histo <- tots <- array(0,dim=c(nx-1,ny-1))
	for(i in 1:(nx-1)){
		for(j in 1:(ny-1)){
		  if(is.null(z)){
			  histo[i,j] <- length(which((x >= xbins[i]) & (x < xbins[i+1]) & (y >= ybins[j]) & (y < ybins[j+1])))
			  tots[i,j] <- length(which((x >= xbins[i]) & (x < xbins[i+1]) & (y >= ybins[j]) & (y < ybins[j+1])))
		  } else{
		    histo[i,j] <- mean(z[which((x >= xbins[i]) & (x < xbins[i+1]) & (y >= ybins[j]) & (y < ybins[j+1]))],na.rm=T)
		    tots[i,j] <- length(which((x >= xbins[i]) & (x < xbins[i+1]) & (y >= ybins[j]) & (y < ybins[j+1])))
		  }
		}
	}
	if(!output_freqs){
		histo <- histo / sum(histo,na.rm=T)
	}
	
	if(output_plot){
		if(exists('image.plot') && imageplot){
			if(!is.null(zlim)){
				image.plot(histo,x=xbins[1:(nx-1)],y=ybins[1:(ny-1)],xlab=xlab,ylab=ylab,axes=axes,zlim=zlim)
			} else{
				image.plot(histo,x=xbins[1:(nx-1)],y=ybins[1:(ny-1)],xlab=xlab,ylab=ylab,axes=axes)
			}
		} else{
			if(!is.null(zlim)){
				image(histo,x=xbins[1:(nx-1)],y=ybins[1:(ny-1)],xlab=xlab,ylab=ylab,zlim=zlim,axes=axes,col=tim.colors(256),bty='o')
			} else{
				image(histo,x=xbins[1:(nx-1)],y=ybins[1:(ny-1)],xlab=xlab,ylab=ylab,axes=axes,col=tim.colors(256),bty='o')
			}
		}
	}
	
	out <- list()
	out$xbins <- xbins
	out$ybins <- ybins
	out$histo <- histo
	out$tots <- tots
	
	invisible(out)
	
}

#FUNCTION: show.distrib.over.time
#Make a color plot showing the distribution of a value over time
#Inputs:
#	values: N x n.times matrix giving the value of a certain measurement for a set of individuals over time (e.g. speed)
#	bins: vector of bins to use for the histogram - defaults to seq(0,1,.05)
#	ts: vector of times - defaults to seq(1,ncol(values))
#Outputs:
#	a color plot visualizing the distribution of values over time
show.distrib.over.time <- function(values,bins=seq(0,1,.05),ts=seq(1,ncol(values))){
	
	#get number of individuals (N) and number of times (n.times)
	N <- nrow(values)
	n.times <- ncol(values)
	
	#create a matrix to hold the distributions for each point in time
	n.bins <- length(bins) #number of bins
	hist.over.time <- matrix(NA,nrow=n.bins,ncol=n.times)
	
	#use a for loop to compute the frequency of values in each bin
	for(i in 1:(n.bins-1)){
		hist.over.time[i,] <- colSums((values >= bins[i]) & (values < bins[i+1]) )
	}
	
	if(exists('image.plot')){
		image.plot(t(hist.over.time),x=ts,y=bins)
	} else{
		image(t(hist.over.time),x=ts,y=bins)
	}
	
}








