#What are the relevant spatial scales associated with animal groups?

#PARAMETERS - MODIFY THESE
general_funcs_path <- '~/Dropbox/code_ari/move_comm_analysis/general/general_function_library.R'
data_path <- '~/Dropbox/hyenas/hyena_data/RData/hyena_xy_level1.RData' #path to xy data
data_path_2 <- '~/Dropbox/hyenas/hyena_data/RData/hyena_timestamps.Rdata' #path to timestamps data, if not included above (else NULL)
R <- 10 #radius (in meters) for spatial headings
speed_dt <- 6 #time difference to use for computing speeds (in units of timesteps, after downsampling if relevant)
downsamp <- 10 #rate of downsampling (to save time), default to 1 to not downsample

#------------------------------------------------
#LIBRARY
library(lubridate)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#FUNCTIONS
#read in functions
source(general_funcs_path)

#LOAD DATA
load(data_path) 
if(!is.null(data_path_2)){
  load(data_path_2)
}

#this if statement is needed because the time vector is called ts in coatis and timestamps in hyenas...
if(exists('ts') & is.POSIXct(ts)){
  timestamps <- ts
}

#DOWNSAMPLE IF NEEDED
if(downsamp > 1){
  n_times_orig <- length(timestamps)
  timestamps <- timestamps[seq(from = 1,to = n_times_orig, by = downsamp)]
  xs <- xs[,seq(from = 1,to = n_times_orig, by = downsamp)]
  ys <- ys[,seq(from = 1,to = n_times_orig, by = downsamp)]
}

#HEADING CORRELATION VS DYADIC DISTANCE

n_inds <- nrow(xs)
n_times <- ncol(xs)
heads <- speeds <- matrix(NA, nrow = n_inds, ncol = n_times)

#break up by day (only really makes sense for discontinuous deployments...)



days <- date(timestamps)
day_start_idxs <- c(1, which(diff(days)==1)+1)

#spatially discretized headings of individuals and speeds
for(i in 1:n_inds){
  for(d in 1:(length(day_start_idxs)-1)){
    tday <- day_start_idxs[d]:(day_start_idxs[d+1]-1)
    headsday <- spatial.headings(xs[i,tday],ys[i,tday],R = R)
    heads[i,tday] <- headsday
    idxs_fut <- tday[seq(speed_dt+1,length(tday),1)]
    idxs_now <- tday[seq(1,length(tday)-speed_dt,1)]
    speeds[i,idxs_now] <- sqrt((xs[i, idxs_fut] - xs[i, idxs_now])^2 + (ys[i, idxs_fut] - ys[i, idxs_now])^2)
  }
}

#dyadic distances and heading correlations
dyad_dists <- dyad_dist_changes <- head_corrs <- speed_diffs <- log_speed_diffs <- array(NA, dim = c(n_inds,n_inds,n_times))
for(i in 1:(n_inds-1)){
  for(j in (i+1):n_inds){
    dyad_dists[i,j,] <- sqrt((xs[i,]-xs[j,])^2 + (ys[i,]-ys[j,])^2)
    head_corrs[i,j,] <- cos(heads[i,])*cos(heads[j,]) + sin(heads[i,])*sin(heads[j,])
    log_speed_diffs[i,j,] <- abs(log(speeds[i,]) - log(speeds[j,]))
    speed_diffs[i,j,] <- abs(speeds[i,] - speeds[j,])
    for(d in 1:(length(day_start_idxs)-1)){
      tday <- day_start_idxs[d]:(day_start_idxs[d+1]-1)
      idxs_fut <- tday[seq(speed_dt+1,length(tday),1)]
      idxs_now <- tday[seq(1,length(tday)-speed_dt,1)]
      dyad_dist_changes[i,j,idxs_now] <- dyad_dists[i,j,idxs_fut] - dyad_dists[i,j,idxs_now]
    }
  }
}

dist_bins <- c(0, 10^seq(0,4,.1))
heads_mean <- speeds_mean <- log_speeds_mean <- change_dyad_dist_mean <- rep(NA, length(dist_bins)-1)
for(i in 1:(length(dist_bins)-1)){
  idxs <- which(dyad_dists >= dist_bins[i] & dyad_dists < dist_bins[i+1])
  heads_mean[i] <- mean(head_corrs[idxs], na.rm=T)
  speeds_mean[i] <- mean(speed_diffs[idxs], na.rm=T)
  log_speeds_mean[i] <- mean(log_speed_diffs[idxs], na.rm=T)
  change_dyad_dist_mean[i] <- mean(dyad_dist_changes[idxs], na.rm=T)
}

#PLOTS

#Plot 1: Distribution of log(dyadic distances) between coatis
quartz()
bins <- seq(-2,3,length=30)
histo <- hist(log(dyad_dists, 10), plot=T, breaks=50, xlab = 'Log (base 10) dyadic distance (m)',main='')

#Plot 2: Mean heading correlation as a function of distance 
quartz()
plot(dist_bins[2:length(dist_bins)],heads_mean, xlab = 'Distance apart (m)', ylab = 'Mean heading correlation', pch = 19, col = '#00000066', cex = 1.5, ylim = c(0,1),log='x')
abline(v=50, col = 'red', lwd = 3, lty = 2)

#Plot 3: Mean difference in speed as a function of distance 
quartz()
plot(dist_bins[2:length(dist_bins)],speeds_mean, xlab = 'Distance apart (m)', ylab = 'Mean absolute speed difference', pch = 19, col = '#00000066', cex = 1.5, ylim = c(0,max(speeds_mean,na.rm=T)), log='x')
abline(h=0, lty = 2)
abline(v=50, col = 'red', lwd = 3, lty = 2)

#Plot 4: Mean change in dyadic distance
quartz()
plot(dist_bins[2:length(dist_bins)],change_dyad_dist_mean, xlab = 'Distance apart (m)', ylab = 'Mean change in dyadic distance per timestep (m)', pch = 19, col = '#00000066', cex = 1.5, ylim = c(min(change_dyad_dist_mean,na.rm=T),max(change_dyad_dist_mean,na.rm=T)), log = 'x')
abline(h=0, lty = 2)
abline(v=50, col = 'red', lwd = 3, lty = 2)
