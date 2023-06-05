#TODO: Fix the border case in event 58, where the lower threshold gets set to 81

#library of general functions

#LIBRARIES
library(dbscan)
library(rgdal)
library(lubridate)
library(stringr)

#GRAPHICS FOR WINDOWS AND MAC
if(.Platform$OS.type=="windows") {
  quartz<-function() windows()
}


#LAT/LON TO UTM CONVERSIONS (AND VICE VERSA)
#Converts a matrix of lons and lats (lons first column, lats second column) to UTM
#Inputs:
#	LonsLats: [N x 2 matrix] of lons (col 1) and lats (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#	EastingsCol1: whether eastings should be given in first column of output (default) or not
#Outputs:
#	EastNorths or NorthEasts: [N x 2 matrix] of Eastings and Northings - eastings are first column by default
latlon.to.utm <- function(LonsLats,EastingsCol1 = TRUE,utm.zone='34',southern_hemisphere=TRUE){
  latlons <- data.frame(X=LonsLats[,2],Y=LonsLats[,1])
  non.na.idxs <- which(!is.na(latlons$X) & !is.na(latlons$Y))
  len <- nrow(latlons)
  non.na.latlons <- latlons[non.na.idxs,]
  coordinates(non.na.latlons) <- ~Y + X
  proj4string(non.na.latlons) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84')
  if(southern_hemisphere){
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +south',sep='')
  } else{
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +north',sep='')
  }
  utm <- spTransform(non.na.latlons,CRS(projection.string))
  EastNorths <- matrix(NA,nrow=len,ncol=2)
  EastNorths[non.na.idxs,] <- utm@coords
  if(!EastingsCol1){
    NorthEasts <- EastNorths[,c(2,1)]
    return(NorthEasts)
  } else{
    return(EastNorths)
  }
}

#Converts a matrix of eastings and northings (eastings first column, northings second column) to latlong
#Inputs:
#	EastNorths: [N x 2 matrix] of eastings (col 1) and northings (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#	LonsCol1: whether lons should be given in first column of output (default) or not
#Outputs:
#	LonLats or LatLons: [N x 2 matrix] of longitudes and latitudes - lons are first column by default 
utm.to.latlon <- function(EastNorths,LonsCol1=TRUE,utm.zone = '34',southern_hemisphere=TRUE){
  utms <- data.frame(X=EastNorths[,1],Y=EastNorths[,2])
  non.na.idxs <- which(!is.na(utms$X) & !is.na(utms$Y))
  len <- nrow(utms)
  if(southern_hemisphere){
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +south',sep='')
  } else{
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +north',sep='')
  }
  non.na.utms <- SpatialPoints(utms[non.na.idxs,],proj4string=CRS(projection.string))
  lonlat <- spTransform(non.na.utms,CRS('+proj=longlat +ellps=WGS84 +datum=WGS84'))
  LonLats <- matrix(NA,nrow=len,ncol=2)
  LonLats[non.na.idxs,] <- lonlat@coords
  if(!LonsCol1){
    LatLons <- LonLats[,c(2,1)]
    return(LatLons)
  } else{
    return(LonLats)
  }	
}


#This function computes information about subgroup membership over time
#Inputs:
# R [numeric]: radius for DBSCAN
# xs, ys [n_inds x n_times matrix]: x and y locations of individuals 
#Outputs:
# subgroup_data: object that contains:
#   $ind_subgroup_membership [n_inds x n_times matrix]: matrix giving the subgroup id for each individual at each time
#   $n_subgroups [n_times vector]: vector of the number of subgroups over time
#   $subgroup_counts [max_n_subgroups x n_times matrix]: matrix giving the number of individuals in each subgroup (with NAs when there are fewer subgroups than the max)
#   $R [numeric]: radius used in DBSCAN
get_subgroup_data <- function(xs, ys, R){
  
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  
  sub_groups <- matrix(NA, nrow = n_inds, ncol = n_times)
  
  for (t_idx in 1:n_times){
    
    #getting the x and y coordinates for each time point
    x <- xs[ , t_idx]
    y <- ys[, t_idx]
    
    non_nas <- which(!is.na(x))
    
    x_non_na <- x[non_nas]
    y_non_na <- y[non_nas]
    
    n_tracked <- length(non_nas)
    
    if(n_tracked >= 1){
      #scan for each time column to get the cluster number for each individual
      groups <- dbscan(x= cbind(x_non_na, y_non_na), eps = R, minPts = 1)
      #saving the cluster for each time into the matrix
      sub_groups[non_nas, t_idx] <- groups$cluster
    }
    
  }
  
  #get max number of groups per time
  ngroups <- suppressWarnings(apply(sub_groups, 2, max, na.rm = T)) #suppress warnings because if all NAs it gives a warning, this is fine because we deal with it in the next step by removing the infinite maxes
  #remove infinite values
  ngroups[is.infinite(ngroups)] <- NA

  
  
  #now want to look at size of the sub-groups as the histogram shows 2 subgroups being most frequent 
  #but this might be because of one individual not being with the group
  
  group_max <- max(ngroups, na.rm = T)
  
  subgroup_counts <- matrix(NA, nrow = group_max, ncol = n_times)
  
  #for loop for getting the number of individuals in each sub group every 10 mins (when the radius for defining within group is 50m, there are maximum 5 groups in full dataset)
  for(i in 1:n_times){
    
    #going through each time stamp
    current_group <- sub_groups[, i]
    #number of individuals in each subgroup 
    counts <- table(current_group)
    #to remove NA's when there is no data from all individuals
    if(length(counts) >= 1){
      #putting the data into the matrix
      subgroup_counts[1:length(counts),i] <- sort(counts)
      
    }
  }
  
  #store all data in a list
  subgroup_data <- list()
  subgroup_data$ind_subgroup_membership <- sub_groups
  subgroup_data$n_subgroups <- ngroups
  subgroup_data$subgroup_counts <- subgroup_counts
  subgroup_data$R <- R
  
  return(subgroup_data)
  
}

#Make a network visualization with colors
#INPUT:
#net[adjacency matrix of the network]
#coati_ids[dataframe of the coati ids]
#OUTPUT:
#plot
visualize_network_matrix <- function(net, coati_ids){
  
  zmin <- min(net, na.rm=T)
  zmax <- max(net, na.rm=T)
  par(mgp=c(3, 1, 0), mar=c(11,11,3,3)) #bottom, left, top, and right ##CHANGE THESE VALUES IF DOING FORLOOP OF THE MATRIX PLOTS AND THEY DON'T FIT mar=c(11,10,6,4))
  image.plot(net, col = viridis(256), zlim=c(zmin,zmax), xaxt= 'n', yaxt = 'n', legend.cex = 7, legend.width = 1.3,legend.mar = 6, axis.args=list(cex.axis=2))
  axis(1, at = seq(0,1,length.out= nrow(net)), labels = coati_ids$name, las = 2, cex.axis=1.8)
  axis(2, at = seq(0,1,length.out= nrow(net)), labels = coati_ids$name, las = 2,  cex.axis=1.8)
  
  points(rep(-.05, nrow(net)),seq(0,1,length.out=n_inds),col=coati_ids$color, xpd = T, pch = 19, cex = 2)
  points(seq(0,1,length.out=nrow(net)),rep(-.05,n_inds),col=coati_ids$color, xpd = T, pch = 19, cex = 2)
}

#changed the size of the labels with cex.axis = 1.5, default is 1
#also changed the size of the legend axis values with: axis.args=list(cex.axis=2), remove if want default
#legend width was 1.3, changed to 5 for the multiple matrix plots
#cex in points change points size
#-0.8 to move the coloured ID points out of the matrix plot
#legend.cex = 5 changed to 7

#for presedente matrix plot, change the points(rep(-.07 for gal to -0.05 for presedente


#------------------------------------------------------------------------

#Function for getting distance between individuals over time
#INPUT:
#xs, ys, radius
#OUTPUT
#object that contains:
  #   $distance over time: distance between each individual at each time point
  #   $proximity_network: probability that i and j are within a distance of R
  #   $r_within [numeric]: radius used for proximity network

get_proximity_data <- function(xs, ys, r_within){

  #get n_inds and n_times
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  
  #get distance between individuals at these subset times
  
  #make array to store data into
  dist_over_time <- array(NA, dim = c(n_inds, n_inds, n_times))
  
  #for loop through each individual with each individual for each time point (where all individuals are together and have gps point)
  for(t in 1:ncol(xs)){
    for(i in 1:n_inds){
      for(j in 1:n_inds){
        
        #get x and y coordinates to calculate distance using pythagorus
        xcoord_i <- xs[i,t]
        xcoord_j <- xs[j,t]
        dx <- (xcoord_i - xcoord_j)
        ycoord_i <- ys[i,t]
        ycoord_j <- ys[j,t]
        dy <- (ycoord_i - ycoord_j)
        dist <- sqrt((dx)^2 + (dy)^2)
        
        #put the statements into the array at the correct position in the correct dimension
        dist_over_time[i, j, t] <- dist
      }
    }
  }
  
  #threshold the distances to determine if individuals are in proximity at each time
  net_over_time <- dist_over_time < r_within
  
  #construct proximity network
  proximity_net <- apply(net_over_time, MARGIN = c(1,2), FUN = mean, na.rm=T)
  diag(proximity_net) <- NA
  
  #store output
  proximity_data <- list()
  proximity_data$proximity_net <- proximity_net
  proximity_data$dist_over_time <- dist_over_time
  proximity_data$r_within <- r_within
  
  #return output
  return(proximity_data)

}


#Randomise the splits and construct new randomized data frame
#Keep original group the same, keep subgroup sizes the same, but randomize who goes to which subgroup
#INPUTS:
# splits_df: data frame with information about splits from real data
#OUTPUTS:
# splits_df_random: randomised version of splits data frame
randomise_splits <- function(splits_df){
  
  #start with the splits_df data frame from the data
  splits_df_rand <- splits_df
  
  #loop over each row and randomize who goes where
  for (i in 1:nrow(splits_df_rand)){
    
    #get the original group and sizes of subgroups for one row
    row <- splits_df[i,]
    orig_group <- row$orig_group[[1]]
    n_sub1 <- row$n_sub1
    n_sub2 <- row$n_sub2
    n_sub3 <- row$n_sub3
    
    #shuffle the original group to a random order
    orig_group_shuffled <- sample(orig_group)
    
    #allocate them to groups
    sub1 <- orig_group_shuffled[1:n_sub1]
    sub2 <- orig_group_shuffled[(length(sub1)+1):(length(sub1)+n_sub2)]
    
    splits_df_rand$sub1[i] <- list(sub1)
    splits_df_rand$sub2[i] <- list(sub2)
    
    if(n_sub3 > 0){
      sub3 <- orig_group_shuffled[(length(sub1)+length(sub2)+1):(n_sub1+n_sub2+n_sub3)]
      splits_df_rand$sub3[i] <- list(sub3)
    }
    
  }
  
  return(splits_df_rand)
  
}

#Helper function for get_consistency
#Returns q if q <= 0.5 and 1-q if q > 0.5
dist_to_0_or_1 <- function(q){
  
  if(is.na(q)){
    return(NA)
  }
  if(q > 0.5) {return(1-q)
   }else{return(q)}
  
}

#Compute our metric of consistency, C
#Let q_dyad = p_dyad_together
#m_dyad = 1-q_dyad if q_dyad > 0.5 or = q_dyad if q_dyad < 0.5 (so it's the distance to 0 or 1, whichever is closer)
#C = 1 - 2*mean(m_dyad)
#At the end, C is a measure of consistency where C = 0 if all q_dyad are at 0.5 and C = 1 if all q_dyad are either 0 or 1
#Why do we multiply by 2 and subtract from 1? The multiplying by 2 puts the number between 0 and 1 instead of between 0 and 0.5. The 1 minus makes it be consistency instead of inconsistency
#INPUTS:
# p_dyad_together: [matrix] of probabilities of splitting together given you were both in the original group
#OUTPUTS:
# C: [numeric] metric of consistency
get_consistency <- function(p_dyad_together){
  
  m_dyad <- sapply(p_dyad_together, FUN = dist_to_0_or_1)
  C <- 1 - 2*mean(m_dyad, na.rm=T)
  return(C)
}

#function that takes in a splits dataframe and outputs a p_dyad_together matrix (probability of splitting together for each dyad)
#INPUTS:
# splits_df_local: [data frame] containing information on group splits
# n_inds_local: [numeric] number of individuals in the full group
#OUTPUTS:
# p_dyad_together_local: [matrix n_inds_local x n_inds_local]: probability of being together in a split, given you were both in the original group 
get_p_dyad_together <- function(splits_df_local, n_inds_local){
  
  #COMPUTE METRIC OF P(STAY TOGETHER | originally together) for each dyad
  p_dyad_together_local <- array(NA, dim = c(n_inds_local,n_inds_local))
  #loop over dyads
  for(i in 1:(n_inds_local-1)){
    for(j in (i+1):n_inds_local){
      
      #find rows where they were both in the original group
      originally_together_rows <- which(unlist(lapply(splits_df_local$orig_group, FUN = function(x){return(i %in% x & j %in% x)})))
      
      #how many times do they end up in the same group
      still_together <- 0
      for(r in originally_together_rows){
        if(i %in% splits_df_local$sub1[r][[1]] & j%in% splits_df_local$sub1[r][[1]]){
          still_together <- still_together + 1
        }
        if(i %in% splits_df_local$sub2[r][[1]] & j %in% splits_df_local$sub2[r][[1]]){
          still_together <- still_together + 1
        }
        if(i %in% splits_df_local$sub3[r][[1]] & j %in% splits_df_local$sub3[r][[1]]){
          still_together <- still_together + 1
        }
        
      }
      
      p_dyad_together_local[i,j] <- still_together / length(originally_together_rows)
      
    }
  }
  return(p_dyad_together_local)
  
}

#Function to match names in fission/fusion manual events table to coati ids idxs
#INPUTS:
# subgroup_names: a character string containing the names of the coatis in a given subgroup
# coati_ids: table of coati ids to match to
#OUTPUTS:
# subgroup_idxs: a vector of the indexes associated with the subgroup members
match_coati_names <- function(subgroup_names, coati_ids){
  names_split <- strsplit(subgroup_names,',')[[1]]
  names_split <- sub(' ','',names_split)
  names_short <- sapply(names_split, function(x){return(substr(x,1,3))})
  
  #match names to get individual indexes
  subgroup_idxs <- match(names_short, coati_ids$name_short)
  return(subgroup_idxs)
  
}

#ANALYSE_FF_EVENT
#Analyse (and, if plot=T make a visualization of) a fission-fusion event
#This function takes in information about a fission or fusion event as well as
#tracking data to identify relevant time points in a fission or fusion event.
#In particular, the start_time and end_time of the event are determined based on
#a double threshold method, then a before_time and after_time are identified. 
#This identifies 3 phases of the event: before (before_time:start_time), during
#(start_time:end_time) and after (end_time:after_time).
#Finally the displacements and speeds of the subgroups during these phases and various 
#relevant angles are calculated. More details below.
#---How are the start and end time identified?---
#As a start, we look at a window of time around an identified event (can be 
#manually or automatically identified). The size of the window is determined by
#the parameter 'max_time', and we go from (tidx - max_time):(tidx + max_time) 
#where tidx is the identified event time index.
#Within the window, we compute the distance between the centroids of the two 
#subgroups that are fissioning or fusioning at each time (again, the subgroups
#could be manually identified or automatically identified). 
#We then use a double threshold method to determine the start and end point of
#the fission or fusion event within that time window. To do so, we first 
#categorize all time points as being above the higher threshold thresh_h (2), 
#between the two thresholds (1), or below the lower threshold thresh_l (0). We 
#then identify contiguous periods of time where the dyadic distance went from 
#2-111111(any number of ones)-0 (i.e. high-mid-low) for a fusion or 0-1111...-2
#(i.e. low-mid-high) for a fission. If there are multiple possible time periods
#detected within the window, we choose the one that is closest to tidx. (See also
#subtlety 1 below).
#---How are the before and after times identified?---
#The before time is defined as the time point tidx - time_window (default 
#time_window = 300s), unless the two groups are not sufficiently together (for
#a fission) or apart (for a fusion) at that time. If the latter, the before_time
#is identified as the point just before the two subgroups cross a threshold midway
#between thresh_l and thresh_h (i.e. at (thresh_l + thresh_h) / 2). The logic here
#is that, for a fusion we are looking for what the full group (combination of the
#two eventual subgroups) was doing before they began to split. However, we do not
#want this point to fall during a prior fusion event, so we require the centroids
#of the two subgroups to be less than (thresh_l + thresh_h) / 2 distance apart.
#The after time is defined in a parallel way, but using the time period after the
#event. It is usually set to tidx + time_window, unless the subgroups come back too
#close together (for a fission) or go too far apart (for a fusion). Again, we use
#the time just prior to crossing the midway point between the two thresholds as
#the end time, if this threshold is crossed before time_window seconds has elapsed.
#---How are the displacements defined?---
#We compute the displacement of the centroid of each subgroup (group A and group B)
#as well as the displacement of the centroid of the combined group (group AB) during
#each of the time windows: before (before_time:start_time), during (start_time:end_time)
#and after (end_time:after_time).
#---How are the angles defined?---
#We define 3 relevant angles relevant to a fission event (might define more later
#for merge events):
# split_angle: this is the angle at which the two groups split. It is defined as
#   the angle traced out by the points p_A(end_time), p_AB(start_time), and p_B(end_time)
#   where p_A(t) is the position of subgroup A's centroid at time (t) and likewise
#   for subgroup B and the combined group. 
# turn_angle_A: the turning angle of subgroup A. This is defined as the angle
#   formed by the points p_AB(before_time), p_AB(start_time), and p_A(end_time)
# turn_angle_B: likewise for subgroup B
# Note that all angles use the point p_AB(start_time) as their central point
#INPUTS:
# events: a data frame of manually-labeled ff events
# i: index to the row to use
# xs, ys: matrices of x and y coordinates [n_inds x n_times]
# ts: vector times in datetime format [length = n_times]
# max_time: [numeric] maximum time forward and back to look for the start and end of the event
# thresh_h: upper threshold for determining when individuals are "apart" (default 50)
# thresh_l: lower threshold for determining when individuals are "together" (default 15)
# time_window: time to move backward or forward in time to identify the before and after times
# plot: [boolean] whether to plot the event or not
#OUTPUTS:
# (if plot == T) a plot showing (top) dyadic distance over time and (bottom) a visualization of trajectories
# out: a list of information extracted about the event
# out$start_time: start time
# out$end_time: end time
# out$before_time: before time
# out$after_time: after time
# out$disps: matrix of displacements of the different subgroups (rows) during the different
#   time intervals (columns). Rows and columns are named for easy access.
# out$speeds: same as disps matrix, but with speeds (so dividing by the time differences)
# out$split_angle: split angle in degrees, description above
# out$turn_angle_A: turn angle for subgroup A in degrees, description above
# out$turn_angle_B: turn angle for subgroup B in degrees, description above.
#SUBTLETIES:
#Subtlety 1: Sometimes, due to the multi-scale nature of these events and the
#fact that we are approximating subgroup locations with centroids, the dyadic
#distance does not go below the lower threshold thresh_l and/or above the upper 
#threshold thresh_h. In this case, we still try to identify the start_time and 
#end_time, but modify the thresholds as follows. 
#First, let's define the period between tidx - max_time and tidx as the prior period,
#the period between tidx  and tidx + max_time as the subsequent period, and the 
#period between (tidx - max_time/2) and (tidx + max_time/2) as the middle period.
#--> For a fission, if the dyadic distance does not go above thresh_h in the subsequent
#period, then we instead replace thresh_h with the maximum - .001 of the dyadic 
#distance during that time period. Second, if the dyadic distance during the middle
#period does not drop below thresh_l, then we instead move thresh_l to the minimum + .001
#of the dyadic distance during that middle period.
#--> For a fusion, everything is reversed. If the dyadic distance does not go
#above thresh_h during the prior period, we replace thresh_h with the maximum - .001
#of the dyadic distance during the subsequent period. And if the dyadic distance does
#not go below thresh_l during the middle period, we replace thresh_l with the minimum + .001
#of the dyadic distance during the middle period.
#In very rare cases, due to these rules the upper bound may get changed to 
#something below the original thresh_l, or the lower bound may get changed to
#something above the original thresh_h. In that case, we revert to the original thresholds.
#Subtlety 2: NA handling. NAs are handled in a couple of different ways.
#--For finding the start_time and end_time, NAs are essentially ignored, and cannot
#be part of the sequence from start_time:end_time (at least not in terms of the)
#dyadic distance... individuals can drop out and they will just be excluded from
#the centroid calculation.
#If no start and end times are found, nothing else is computed (all vaules are
#filled in with NAs or NULL for matrices).
#--For finding the before_time and after_time, if an NA is hit in the forward or
#backward direction, the before_time (or respectively, the after_time) is marked
#as NA. Metrics involving that time point can then not be computed and are also 
#NA. 
#--If the start and end time are on different days, then both are given NA.
#--If the before and start time, or after and end time, are on different days, then
#both are given NA. 
#--If the before or after time are NA, then other metrics stemming from those times get NA
analyse_ff_event <- function(i, events, xs, ys, ts, max_time = 600, thresh_h = 50, thresh_l = 15, time_window = 300, plot = T){
  t_event <- events$tidx[i] #time of the event
  group_A <- events$group_A_idxs[i][[1]] #group A individual idxs
  group_B <- events$group_B_idxs[i][[1]] #group B individual idxs
  group_A_names <- events$group_A[i]
  group_B_names <- events$group_B[i]
  ti <- t_event - max_time #initial time to plot
  tf <- t_event + max_time #final time to plot
  
  if(ti <1 | tf > length(ts)){
    out <- list()
    out$start_time <- out$end_time <- out$before_time <- out$after_time <- NA
    out$disps <- out$speeds <- NULL
    out$turn_angle_A <- out$turn_angle_B <- out$split_angle <- NA
    return(out)
  }
  
  event_type <- events$event_type[i]
  datetime <- events$datetime[i]
  nA <- length(group_A)
  nB <- length(group_B)
  nT <- ncol(xs)
  
  #get x and y coordinates of the relevant individuals in each subgroup
  xA <- matrix(xs[group_A,],nrow=nA,ncol=nT)
  xB <- matrix(xs[group_B,],nrow=nB,ncol=nT)
  yA <- matrix(ys[group_A,],nrow=nA,ncol=nT)
  yB <- matrix(ys[group_B,],nrow=nB,ncol=nT)
  
  #get centroids of the two subgroups
  xcA <- colMeans(xA, na.rm=T)
  ycA <- colMeans(yA, na.rm=T)
  xcB <- colMeans(xB, na.rm=T)
  ycB <- colMeans(yB, na.rm=T)
  
  #get distance between centroids
  dyad_dist <- sqrt((xcA - xcB)^2 + (ycA - ycB)^2)
  
  #classify the dyadic distance into categories:
  #0 = below lower thershold
  #1 = between thresholds
  #2 = above higher threshold
  dyad_dist_event <- dyad_dist[ti:tf]
  
  #first consider modifying thresholds according to subtlety 1 above
  upper <- thresh_h
  lower <- thresh_l
  after_idxs <- (max_time+1):(2*max_time+1) #indexes after the marked event
  middle_idxs <- (max_time / 2):(max_time*3/2)
  before_idxs <- 1:max_time #indexes before the marked event
  if(event_type == 'fission'){
    if(sum(!is.na(dyad_dist_event[after_idxs]))>1){
      if(max(dyad_dist_event[after_idxs],na.rm=T) < thresh_h){
        upper <- max(dyad_dist_event[after_idxs],na.rm=T) - .001
      } else{
        upper <- thresh_h
      }
    } 
    if(sum(!is.na(dyad_dist_event[middle_idxs]))){
      if(min(dyad_dist_event[middle_idxs],na.rm=T) > thresh_l){
        lower <- min(dyad_dist_event[middle_idxs],na.rm=T) + .001
      } else{
        lower <- thresh_l
      }
    }
  }
  if(event_type == 'fusion'){
    if(sum(!is.na(dyad_dist_event[before_idxs]))>1){
      if(max(dyad_dist_event[before_idxs],na.rm=T) < thresh_h){
        upper <- max(dyad_dist_event[before_idxs],na.rm=T) - .001
      } else{
        upper <- thresh_h
      }
    }
    if(sum(!is.na(dyad_dist_event[middle_idxs]))>1){
      if(min(dyad_dist_event[middle_idxs],na.rm=T) > thresh_l){
        lower <- min(dyad_dist_event[middle_idxs],na.rm=T) + .001
      } else{
        lower <- thresh_l
      }
    }
  }
  
  #if the upper bound was changed to something < thresh_l, move threshold back to thresh_h
  if(upper <= thresh_l){
    upper <- thresh_h
  }
  #likewise for lower bound
  if(lower >= thresh_h){
    lower <- thresh_l
  }
  
  #get category of each moment in time
  #0 = below lower, 1 = middle, 2 = above upper
  category <- rep(NA, length(dyad_dist_event))
  category[which(dyad_dist_event < lower)] <- 0
  category[which(dyad_dist_event >= lower & dyad_dist_event < upper)] <- 1
  category[which(dyad_dist_event >= upper)] <- 2
  category[which(is.na(dyad_dist_event))] <- 3 #NAs are denoted with 3
  
  #run length encoding to get sequences of each category
  seqs <- rle(category)
  
  #find sequences of high-middle-low (2,1,0) or low-mid-high (0,1,2)
  seqs_str <- paste0(as.character(seqs$values), collapse = '') #convert to string
  if(event_type=='fission'){
    event_loc <- as.data.frame(str_locate_all(seqs_str,'012')[[1]])
  } 
  if(event_type == 'fusion'){
    event_loc <- as.data.frame(str_locate_all(seqs_str,'210')[[1]])
  }
  
  #for seuqneces of hml or lmh (for fission and fusion respectively), get the time index when they start and end
  if(nrow(event_loc)>0){
    for(r in 1:nrow(event_loc)){
      event_loc$start_time[r] <- ti + sum(seqs$lengths[1:event_loc$start[r]])
      event_loc$end_time[r] <- ti + sum(seqs$lengths[1:event_loc$end[r]-1])
    }
  }
 
  if(plot == T){
    quartz() #open a new plot for mac
    par(mfrow=c(2,1))
    plot(ti:tf, dyad_dist[ti:tf],type='l', main = paste(event_type, datetime),xlab='Time (min)',ylab = 'Distance apart (m)')
    abline(v=t_event,col='black', lty = 2)
    abline(h = thresh_h, col = 'darkorange1')
    abline(h = thresh_l, col = 'magenta')
    for(r in 1:nrow(event_loc)){
      abline(v=event_loc$start_time[r], col = 'green')
      abline(v=event_loc$end_time[r], col = 'green')
    }
    
    xmin <- min(min(xA[,ti:tf],na.rm=T),min(xB[,ti:tf],na.rm=T))
    xmax <- max(max(xA[,ti:tf],na.rm=T),max(xB[,ti:tf],na.rm=T))
    ymin <- min(min(yA[,ti:tf],na.rm=T),min(yB[,ti:tf],na.rm=T))
    ymax <- max(max(yA[,ti:tf],na.rm=T),max(yB[,ti:tf],na.rm=T))
    plot(NULL, xlim=c(xmin,xmax),ylim=c(ymin,ymax),asp=1, xlab='Easting', ylab = 'Northing', main = paste('(Red =', group_A_names, '), (Blue =', group_B_names,')'))
    for(j in 1:nrow(xA)){
      lines(xA[j,ti:tf],yA[j,ti:tf],type='l',col='#FF000033')
    }
    for(j in 1:nrow(xB)){
      lines(xB[j,ti:tf],yB[j,ti:tf],type='l',col='#0000FF33')
    }
    lines(xcA[ti:tf],ycA[ti:tf], col = '#FF0000', lwd = 2)
    lines(xcB[ti:tf],ycB[ti:tf], col = '#0000FF', lwd = 2)
    
    points(xcA[t_event], ycA[t_event], pch = 8, col = 'black')
    points(xcB[t_event], ycB[t_event], pch = 8, col = 'black')
    points(xcA[ti], ycA[ti], pch = 1, col = 'black')
    points(xcB[ti], ycB[ti], pch = 1, col = 'black')
    points(xcA[tf], ycA[tf], pch = 4, col = 'black')
    points(xcB[tf], ycB[tf], pch = 4, col = 'black')
    
    #algorithm-identified start and end
    points(xcA[event_loc$start_time],ycA[event_loc$start_time],pch = 1, col = 'green')
    points(xcB[event_loc$start_time],ycB[event_loc$start_time],pch = 1, col = 'green')
    points(xcA[event_loc$end_time],ycA[event_loc$end_time],pch = 4, col = 'green')
    points(xcB[event_loc$end_time],ycB[event_loc$end_time],pch = 4, col = 'green')
    
  }
  
  #create an object to output the start and end times
  start_time <- event_loc$start_time
  end_time <- event_loc$end_time
  
  #if there is more than one start time, go with the closest to the identified fission or fusion point
  if(length(start_time)>1){
    ff_time <- events$tidx[i]
    time_diff <- abs(start_time-ff_time)
    start_time <- start_time[which(time_diff==min(time_diff))]
  }
  #same for end time
  if(length(end_time)>1){
    ff_time <- events$tidx[i]
    time_diff <- abs(end_time-ff_time)
    end_time <- end_time[which(time_diff==min(time_diff))]
  }
  
  #if start time not found (NULL) change to NA
  if(is.null(start_time)){start_time <- NA}
  if(is.null(end_time)){end_time <- NA}
  
  #if the start time is on a different date from the end time, make both NA
  if(!is.na(start_time) & !is.na(end_time)){
    if(date(ts[start_time])!= date(ts[end_time])){
      start_time <- NA
      end_time <- NA
    }
  }
  
  #save to output
  out <- list(start_time = start_time, end_time = end_time)
  
  if(!is.na(start_time) & !is.na(end_time)){
    #GET BEFORE AND AFTER TIMES
    #find the before_time and after_time (times before the start time and after 
    #the end time of the event)
    min_before_time <- start_time - time_window
    max_after_time <- end_time + time_window
    thresh_m <- (thresh_h + thresh_l)/2 #middle threshold is average of upper and lower
    #go backward in time until the two groups cross thresh_m or until the time window has elapsed 
    for(t in seq(start_time, min_before_time, -1)){
      xc_A <- mean(xs[events$group_A_idxs[i][[1]], t], na.rm=T)
      yc_A <- mean(ys[events$group_A_idxs[i][[1]], t], na.rm=T)
      xc_B <- mean(xs[events$group_B_idxs[i][[1]], t], na.rm=T)
      yc_B <- mean(ys[events$group_B_idxs[i][[1]], t], na.rm=T)
      dist_apart <- sqrt((xc_A - xc_B)^2 + (yc_A - yc_B)^2)
      #if you hit an NA, then the before time is undefined - t <- NA and break
      if(is.na(dist_apart)){
        t <- NA
        break
      }
      if(events$event_type[i] == 'fission' & dist_apart > thresh_m){
        break
      }
      if(events$event_type[i] == 'fusion' & dist_apart < thresh_m){
        break
      }
    }
    #store the time, which will either be determined by the time window or by the 
    #groups crossing the threshold
    before_time <- t
    
    #if the before time is on a different date, make it NA
    if(!is.na(before_time)){
      if(date(ts[before_time])!= date(ts[start_time])){
        before_time <- NA
      }
    }
    
    #same thing for the after times
    for(t in seq(end_time, max_after_time, 1)){
      xc_A <- mean(xs[events$group_A_idxs[i][[1]], t], na.rm=T)
      yc_A <- mean(ys[events$group_A_idxs[i][[1]], t], na.rm=T)
      xc_B <- mean(xs[events$group_B_idxs[i][[1]], t], na.rm=T)
      yc_B <- mean(ys[events$group_B_idxs[i][[1]], t], na.rm=T)
      dist_apart <- sqrt((xc_A - xc_B)^2 + (yc_A - yc_B)^2)
      #if you hit an NA, then the after time is undefined - t <- NA and break
      if(is.na(dist_apart)){
        t <- NA
        break
      }
      if(events$event_type[i] == 'fission' & dist_apart < thresh_m){
        break
      }
      if(events$event_type[i] == 'fusion' & dist_apart > thresh_m){
        break
      }
    }
    
    after_time <- t
    
    #if the after time is on a different date, make it NA
    if(!is.na(after_time)){
      if(date(ts[after_time])!= date(ts[end_time])){
        after_time <- NA
      }
    }
    
    #store the before and after times for later returning
    out$before_time <- before_time
    out$after_time <- after_time
    
    #GET SPEEDS BEFORE AND AFTER 
    #Get centroid locations for time_before, time_start, time_end, time_after
    times <- c(before_time, start_time, end_time, after_time)
    xs_AB <- ys_AB <- matrix(NA, nrow = 3, ncol = length(times))
    row.names(xs_AB) <- row.names(ys_AB) <- c('A','B','AB')
    colnames(xs_AB) <- colnames(ys_AB) <- c('before_time','start_time','end_time','after_time')
    for(t in 1:length(times)){
      if(!is.na(times[t])){
        xs_AB[1,t] <- mean(xs[events$group_A_idxs[i][[1]], times[t]],na.rm=T)
        xs_AB[2,t] <- mean(xs[events$group_B_idxs[i][[1]], times[t]],na.rm=T)
        xs_AB[3,t] <- mean(xs[c(events$group_A_idxs[i][[1]],events$group_B_idxs[i][[1]]), times[t]],na.rm=T)
        ys_AB[1,t] <- mean(ys[events$group_A_idxs[i][[1]], times[t]],na.rm=T)
        ys_AB[2,t] <- mean(ys[events$group_B_idxs[i][[1]], times[t]],na.rm=T)
        ys_AB[3,t] <- mean(ys[c(events$group_A_idxs[i][[1]],events$group_B_idxs[i][[1]]), times[t]],na.rm=T)
      }
    }
    
    #get displacements before during and after
    xdiffs <- t(diff(t(xs_AB)))
    ydiffs <- t(diff(t(ys_AB)))
    
    colnames(xdiffs) <- colnames(ydiffs) <- c('before','during','after')
    disps <- sqrt(xdiffs^2 + ydiffs^2)
    out$disps <- disps
    
    #speeds in m / min
    dts <- diff(times)
    speeds <- disps
    for(j in 1:length(dts)){
      speeds[,j] <- disps[,j] / dts[j] * 60
    }
    
    out$speeds <- speeds
    
    #GET ANGLES
    #split angle, turn angle of A, turn angle of B relative to initial group heading
    #see diagram
    split_angle <- angle_between_vectors(x1_i = xs_AB['AB','start_time'],
                                         y1_i = ys_AB['AB','start_time'],
                                         x1_f = xs_AB['A','end_time'],
                                         y1_f = ys_AB['A', 'end_time'],
                                         x2_i = xs_AB['AB', 'start_time'],
                                         y2_i = ys_AB['AB', 'start_time'],
                                         x2_f = xs_AB['B', 'end_time'],
                                         y2_f = ys_AB['B', 'end_time'])
    turn_angle_A <- angle_between_vectors(x1_i = xs_AB['AB','before_time'],
                                          y1_i = ys_AB['AB','before_time'],
                                          x1_f = xs_AB['AB','start_time'],
                                          y1_f = ys_AB['AB', 'start_time'],
                                          x2_i = xs_AB['AB', 'start_time'],
                                          y2_i = ys_AB['AB', 'start_time'],
                                          x2_f = xs_AB['A', 'end_time'],
                                          y2_f = ys_AB['A', 'end_time'])
    turn_angle_B <- angle_between_vectors(x1_i = xs_AB['AB','before_time'],
                                          y1_i = ys_AB['AB','before_time'],
                                          x1_f = xs_AB['AB','start_time'],
                                          y1_f = ys_AB['AB', 'start_time'],
                                          x2_i = xs_AB['AB', 'start_time'],
                                          y2_i = ys_AB['AB', 'start_time'],
                                          x2_f = xs_AB['B', 'end_time'],
                                          y2_f = ys_AB['B', 'end_time'])
    
    out$turn_angle_A <- turn_angle_A
    out$turn_angle_B <- turn_angle_B
    out$split_angle <- split_angle
  } else{
    out$start_time <- out$end_time <- out$before_time <- out$after_time <- NA
    out$turn_angle_A <- out$turn_angle_B <- out$split_angle <- NA
    out$disps <- out$speeds <- NULL
  }
  
  #return the output object
  invisible(out)
}

#Angle between vectors in degrees
angle_between_vectors <- function(x1_i, y1_i, x1_f, y1_f, x2_i, y2_i, x2_f, y2_f){
  
  dx1 <- x1_f - x1_i
  dx2 <- x2_f - x2_i
  dy1 <- y1_f - y1_i
  dy2 <- y2_f - y2_i
  
  #get dot product
  dot <- dx1*dx2 + dy1*dy2
  
  #magnitude of vectors (length)
  mag1 <- sqrt(dx1^2 + dy1^2)
  mag2 <- sqrt(dx2^2 + dy2^2)
  
  #cos of angle
  cosang <- dot / (mag1* mag2)
  
  #get the angle
  angle <- acos(cosang)*180/pi
  
  return(angle)

}

#DETECT FF EVENTS
#Detect fissions and fusions using "sticky-DBSCAN" method from Libera et al. 
#Start by defining an adjacency matrix ('together' in the code) of which dyads 
#are "connected" at any moment in time. Dyads are considered to be connected if
#they go within a distance R_inner of one another, and they continue to be 
#connected until they leave a distance R_outer of one another on both ends
#(before and after) of the period where their distance dropped below R_inner. 
#This double threshold makes periods of connectedness more stable by removing 
#the "flicker" that would result from having a single threshold.
#NAs are handled by ending a period of connectedness if an individual has an NA
#at the point immediately before / after (if the were connected after / before 
#that NA). Individuals with NAs are not included in the together matrix and will
#not be included in the groups.
#Once connectedness of dyads is determined, merge dyads together into groups by
#using DBSCAN on 1 - together as the distance matrix, with eps = something small
#(.01 in the code). Store these groups in groups_list, a list of lists whose 
#length is equal to n_times.
#Stepping through the groups_list, identify all changes in group membership,
#i.e. consecutive time points when the groups do not match. The algorithm 
#flags any change, including instances where individuals disappear or reappear 
#in groups due to missing data (but these are later ignored). Store in "changes"
#data frame.
#In a last step, find all changes where the individuals were preserved (i.e. 
#individuals did not disappear) and where the number of subgroups either 
#increased (a fission) or decreased (a fusion). Store these in a data frame 
#called events_detected.
#In cases where the number of subgroups went from 1 to 2 (for a fission) or 2 to
#1 (for a fusion), identify the members of the 2 subgroups as group A and B 
#(labels arbitrary). In cases where there were more subgroups, first 
#INPUTS:
# R_inner, R_outer: [numeric] inner and outer thresholds respectively
# xs, ys: [n_inds x n_times matrices] of x and y UTM coordinates
# ts: vector of timestamps
# coati_ids: coati ids data frame
#OUTPUTS:
# out: a list of outputs containing
#   out$events_detected: data frame with info on detected fissions and fusions.
#     $tidx: (initial) time index of the event
#     $event_type: "fissin" or "fusion"
#     $group_A_idxs, $group_B_idxs: individual idxs of subgroup members
#     $group_A, $group_B: first 3 letters of names of subgroup members
#     $n_A, n_B: number of individuals in each subgroup
#   out$groups_list: list of subgroups in each timestep
#     groups_list[[t]] gives a list of the subgroups
#     groups_list[[t]][[1]] gives the vector of the first subgroup, etc.
#   out$together: [n_inds x n_inds x n_times matrix] of whether individuals are 
#     connected (1) or not (0) or unknown (NA)
#   changes: data frame containing all the subgroups membership changes (not 
#     just ones identified as fissions or fusions)
detect_fissions_and_fusions <- function(R_inner, R_outer, xs = xs, ys = ys, ts = ts, coati_ids = coati_ids, verbose = T){
  
  #----Identify subgroups at each point
  if(verbose){print('Identifying subgroups at each point using sticky DBSCAN')}
  #number of inds and times
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  
  #day start indexes
  days <- date(ts)
  day_start_idxs <- c(1, which(diff(days)==1)+1)
  
  #Get dyadic distances for each pair, then use double threhsold method to determine if htey are otgether at any moment
  dyad_dists <- together <- array(NA, dim = c(n_inds, n_inds, n_times))
  for(i in 1:(n_inds-1)){
    for(j in (i+1):n_inds){
      
      #dyadic distance
      dx <- xs[i,] - xs[j,]
      dy <- ys[i,] - ys[j,]
      dyad_dists[i,j,] <- sqrt(dx^2 + dy^2)
      
      #together or not
      #loop over days. for each day...
      for(d in 1:(length(day_start_idxs)-1)){
        
        #get times for that day
        t_day <- day_start_idxs[d]:(day_start_idxs[d+1]-1)
        
        #dyadic distances on that day
        dyad_dists_ij <- dyad_dists[i,j,t_day]
        
        if(sum(!is.na(dyad_dists_ij))==0){
          next
        }
        
        #times when together within inner radius
        together_inner <- dyad_dists_ij <= R_inner
        
        #times when together within outer radius
        together_outer <- dyad_dists_ij <= R_outer
        
        together_ij <- together_inner
        
        if(sum(together_inner,na.rm=T)==0){
          together[i,j,t_day] <- together[j,i,t_day] <- together_ij
          next
        }
        
        #go backwards from crossing points into inner radius to find the 'starts' when crossed the outer radius
        inner_starts <- which(diff(together_inner)==1)+1
        if(length(inner_starts)==0){
          together[i,j,t_day] <- together[j,i,t_day] <- together_ij
          next
        }
        for(k in 1:length(inner_starts)){
          crossing <- inner_starts[k]
          curr_time <- crossing
          for(curr_time in seq(crossing,1,-1)){
            if(is.na(together_outer[curr_time])){
              start <- curr_time + 1
              break
            }
            if(curr_time == 1){
              start <- curr_time
              break
            }
            if(together_outer[curr_time]==F){
              start <- curr_time + 1
              break
            }
            
          }
          together_ij[start:crossing] <- T
        }
        
        #go forwards from crossing points out of outer radius to find the 'ends' when crossed the outer radius
        inner_ends <- which(diff(together_inner)==-1)
        if(length(inner_ends)==0){
          together[i,j,t_day] <- together[j,i,t_day] <- together_ij
          next
        }
        for(k in 1:length(inner_ends)){
          crossing <- inner_ends[k]
          curr_time <- crossing
          for(curr_time in seq(crossing,length(together_ij),1)){
            if(is.na(together_outer[curr_time])){
              end <- curr_time - 1
              break
            }
            if(curr_time == length(together_outer)){
              end <- curr_time
              break
            }
            if(together_outer[curr_time]==F){
              end <- curr_time - 1
              break
            }
            
          }
          together_ij[crossing:end] <- T
        }
        
        together[i,j,t_day] <- together[j,i,t_day] <- together_ij
        
      }
    }
  }
  
  #Identify groups from together matrices
  groups <- matrix(NA, nrow = n_inds, ncol = n_times)
  for(t in 1:n_times){
    non.nas <- which(colSums(!is.na(together[,,t]))>0)
    if(length(non.nas)<=1){
      next
    }
    non.nas.together <- together[non.nas,non.nas,t]
    diag(non.nas.together) <- 1
    grps.non.nas <- dbscan(x = as.dist(1 - non.nas.together), eps = .1,minPts=1)$cluster
    groups[non.nas, t] <- grps.non.nas
  }
  
  #store groups as lists of lists
  groups_list <- list()
  for(t in 1:n_times){
    
    #create a list to hold the groups at that timestep
    groups_list[[t]] <- list()
    if(sum(!is.na(groups[,t]))==0){
      next
    }
    
    #
    max_group_id <- max(groups[,t],na.rm=T)
    for(i in seq(1,max_group_id)){
      groups_list[[t]][[i]] <- which(groups[,t]==i)
    }
    
  }
  
  #Identifying changes in group membership in consecutive time steps
  
  if(verbose){print('Identifying changes in group membership')}
  event_times <- c()
  for(d in 1:(length(day_start_idxs)-1)){
    t_day <- day_start_idxs[d]:(day_start_idxs[d+1]-1)
    for(t in t_day[1:length(t_day)]){
      groups_curr <- groups_list[[t]]
      groups_next <- groups_list[[t+1]]
      if(sum(!is.na(groups_curr))==0 | sum(!is.na(groups_next))==0){
        next
      }
      n_groups_curr <- length(groups_curr)
      n_groups_next <- length(groups_next)
      matches <- rep(F, n_groups_curr)
      for(i in 1:n_groups_curr){
        group_curr <- groups_curr[[i]]
        matched <- F
        for(j in 1:n_groups_next){
          group_next <- groups_next[[j]]
          if(setequal(group_curr,group_next)){
            matched <- T
          }
        }
        matches[i] <- matched
      }
      if(sum(matches)<n_groups_curr){
        event_times <- c(event_times, t)
      }
    }
  }
  
  if(verbose){print('Finding fissions and fusion')}
  #go through event times and classify events into types
  changes <- data.frame(tidx = event_times)
  changes$datetime <- ts[changes$tidx]
  changes$event_type <- NA
  changes$n_groups_curr <- unlist(lapply(groups_list, length))[event_times]
  changes$n_groups_next <- unlist(lapply(groups_list, length))[event_times+1]
  changes$n_inds_curr <- changes$n_inds_next <- NA
  for(i in 1:nrow(changes)){
    t <- changes$tidx[i]
    groups_curr <- groups_list[[t]]
    groups_next <- groups_list[[t+1]]
    changes$n_inds_curr[i] <- sum(unlist(lapply(groups_curr, length)))
    changes$n_inds_next[i] <- sum(unlist(lapply(groups_next, length)))
  }
  
  changes$event_type[which(changes$n_groups_curr < changes$n_groups_next & changes$n_inds_curr==changes$n_inds_next)] <- 'fission'
  changes$event_type[which(changes$n_groups_curr > changes$n_groups_next & changes$n_inds_curr==changes$n_inds_next)] <- 'fusion'
  
  #get groups for fissions and fusions
  changes$group_A_idxs <- changes$group_A <- changes$group_B_idxs <- changes$group_B <- list(c(0))
  fissions <- which(changes$event_type=='fission')
  fusions <- which(changes$event_type=='fusion')
  
  #Identify which subgroups split or merged
  if(verbose){print('Determinig which subgroups split and merged')}
  for(i in fissions){
    
    groups_curr <- groups_list[[changes$tidx[i]]]
    groups_next <- groups_list[[changes$tidx[i]+1]]
    n_groups_next <- changes$n_groups_next[i]
    n_groups_curr <- changes$n_groups_curr[i]
    
    #if there are only 2 groups at the end, then name them A and B
    if(n_groups_next==2){
      changes$group_A_idxs[i] <- groups_next[1]
      changes$group_B_idxs[i] <- groups_next[2]
    }
    
    #if there are more than 2 groups, need to figure out which one changed
    if(n_groups_next > 2){
      matched_groups <- c()
      for(j in 1:n_groups_next){
        matched <- F
        for(k in 1:n_groups_curr){
          if(setequal(groups_next[[j]], groups_curr[[k]])){
            matched <- T
          }
        }
        if(matched){
          matched_groups <- c(matched_groups, j)
        }
      }
      unmatched_groups <- setdiff(1:n_groups_next, matched_groups)
      if(length(unmatched_groups)==2){
        changes$group_A_idxs[i] <- groups_next[unmatched_groups[1]]
        changes$group_B_idxs[i] <- groups_next[unmatched_groups[2]]
      }else{
        warning(paste('Could not identify unambiguously subgroups for event',i,
                      'due to more than 2 unmatched subgroups'))
      }
      
    }
  }
  
  for(i in fusions){
    
    groups_curr <- groups_list[[changes$tidx[i]]]
    groups_next <- groups_list[[changes$tidx[i]+1]]
    n_groups_next <- changes$n_groups_next[i]
    n_groups_curr <- changes$n_groups_curr[i]
    
    #if there are only 2 groups at the start, then name them A and B
    if(n_groups_curr==2){
      changes$group_A_idxs[i] <- groups_curr[1]
      changes$group_B_idxs[i] <- groups_curr[2]
    }
    
    #if there are more than 2 groups, need to figure out which one changed
    if(n_groups_curr > 2){
      matched_groups <- c()
      for(j in 1:n_groups_curr){
        matched <- F
        for(k in 1:n_groups_next){
          if(setequal(groups_curr[[j]], groups_next[[k]])){
            matched <- T
          }
        }
        if(matched){
          matched_groups <- c(matched_groups, j)
        }
      }
      unmatched_groups <- setdiff(1:n_groups_curr, matched_groups)
      if(length(unmatched_groups)==2){
        changes$group_A_idxs[i] <- groups_curr[unmatched_groups[1]]
        changes$group_B_idxs[i] <- groups_curr[unmatched_groups[2]]
      }else{
        warning(paste('Could not identify unambiguously subgroups for event',i,
                      'due to more than 2 unmatched subgroups'))
      }
      
    }
  }
  
  #get coati names
  coati_ids$name_short <- sapply(coati_ids$name, function(x){return(substr(x,1,3))})
  for(i in c(fissions,fusions)){
    changes$group_A[i] <- list(coati_ids$name_short[changes$group_A_idxs[i][[1]]])
    changes$group_B[i] <- list(coati_ids$name_short[changes$group_B_idxs[i][[1]]])
  }
  
  #final events dataframe (call it events_detected)
  events_detected <- changes[c(fissions, fusions),c('tidx','datetime','event_type','group_A_idxs','group_B_idxs','group_A','group_B')]
  
  #get number of individuals in each subgroup
  events_detected$n_A <- sapply(events_detected$group_A_idxs, length)
  events_detected$n_B <- sapply(events_detected$group_B_idxs, length)
  
  #sort by time
  events_detected <- events_detected[order(events_detected$tidx),]
  
  #add index column
  events_detected$event_idx <- 1:nrow(events_detected)
  
  if(verbose){print('Done.')}
  
  #return things
  out <- list(events_detected = events_detected, groups_list = groups_list, together = together, changes = changes, R_inner = R_inner, R_outer = R_outer)
  return(out) 
  
}



#Spatially discretized headings function from move_comm_analysis/audio_gps_processing/spatially_discretized_headings.R in Ari's github

#Function to get 'spatially discretized heading' for a given trajectory (x,y) 
#This is defined as the vector point from the individual's current location to
#its location after it mas moved a distance of at least R.
#INPUTS:
# x: vector of x coordinates for the trajectory
# y: vector of y coordinates for the trajectory
# R: radius used to compute the headings
# t.idxs: time indexes at which to compute the headings (defaults to entire trajectory)
# backward: boolean, represents whether to go backward in time from current position (if T), or forward (if F, default)
#OUTPUTS:
# spat.heads: vector of spatially discretized headings
spatial.headings <- function(x,y,R,t.idxs=1:length(x),backward=F){
  
  #initialize
  n.times <- length(x)
  spat.heads <- rep(NA,n.times)
  
  #go backwards for backwards vectors
  if(backward){
    t.idxs <- rev(t.idxs)
  }
  
  #loop over all times
  for(t in t.idxs){
    
    #get current location
    x0 <- x[t]
    y0 <- y[t]
    
    if(is.na(x0)){
      next
    }
    
    #move forward (or backward) until radius reached
    found <- 0
    na.found <- 0
    time.vec <- t:n.times
    if(backward){
      time.vec <- seq(t,1,-1)
    }
    for(i in time.vec){
      
      if(backward){
        dx <- x0 - x[i]
        dy <- y0 - y[i]
      } else{
        dx <- x[i] - x0
        dy <- y[i] - y0
      }
      dist <- sqrt(dx^2+dy^2)
      
      if(is.na(dist)){
        spat.heads[t] <- NA
        found <- 1
        na.found <- 1
        break
      }
      else{
        if(dist >= R){
          found <- 1
          break
        }
      }
    }
    
    #if you reach the end of the trajectory, return (and leave rest as NAs)
    if(found){
      if(!na.found){
        spat.heads[t] <- atan2(dy,dx)
      }
    } else{
      return(spat.heads)
    }
  }
  return(spat.heads)
}





 