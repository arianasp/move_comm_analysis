#What are the relevant spatial scales associated with dyadic interactions?

#This script produces the following plots for each species:

#1. Histogram of log dyadic distance between pairs of hyenas (relevant spatial scales)
#2. Mean heading correlation as a function of dyadic distance between hyenas (directional coordination)
#3. Mean difference in speed between pairs of hyenas as a function of distance between them (speed coordination)
#4. Mean change in dyadic distance (over the course of speed_dt seconds) as a function of original dyadic distance

#PARAMS
heading_subsamp <- 60 #subsample rate for headings - to speed up computation
heading_R <- 10 #spatial discretization value for heading computations
speed_dt <- 60 #time step for computing individual speed
fpt_thresh <- 600 #threshold of first passage time over which a heading is not computed (because the individual is not moving so has undefined heading)
largest_bin_quantile <- 0.9999 #99.99% quantile distance apart is the largest bin
begin_day_remove_sec <- 600 #number of seconds to remove at the beginning of each day (to avoid day breaks)
n_bins <- 30 #number of dist apart bins

#DATA
#paths to the data files to use
short_names <- c('Meerkat_HM2017',
                 'Meerkat_HM2019',
                 'Meerkat_L2019',
                 'Meerkat_RW2021',
                 'Meerkat_NQ2021',
                 'Meerkat_ZU2021',
                 'Meerkat_SI2022',
                 'Coati_GA2022',
                 'Coati_PR2022', 
                 'Hyena_TW2017',
                 'Hyena_SS2023'
                 )
files <- c('~/Dropbox/meerkats/processed_data/vocal_interactions_paper_data_submitted/HM2017_COORDINATES_all_sessions.RData',
           '~/Dropbox/meerkats/processed_data/vocal_interactions_paper_data_submitted/HM2019_COORDINATES_all_sessions.RData',
           '~/Dropbox/meerkats/processed_data/vocal_interactions_paper_data_submitted/L2019_COORDINATES_all_sessions.RData',
           '~/Dropbox/meerkats/processed_data/RW2021_COORDINATES_all_sessions.RData',
           '~/Dropbox/meerkats/processed_data/NQ2021_COORDINATES_all_sessions.RData',
           '~/Dropbox/meerkats/processed_data/ZU2021_COORDINATES_all_sessions.RData',
           '~/Dropbox/meerkats/processed_data/SI2022_COORDINATES_all_sessions.RData',
           '~/Dropbox/coati/processed/galaxy/galaxy_xy_highres_level1.RData',
           '~/Dropbox/coati/processed/presedente/presedente_xy_highres_level1.RData',
           '~/Dropbox/hyenas/hyena_data/RData/hyena_xy_level1.RData',
           '~/Dropbox/Tag-A-Clan2023/processed/hyena_south2023_xy_level0.RData'
           )

#LIBRARY
library(lubridate)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#DIRECTORIES AND PARAMETERS
codedir <- '~/Dropbox/code_ari/hyena_fission_fusion/'
plotdir <- '~/Dropbox/cross_species/spatial_scales/'

#FUNCTIONS
#read in functions
setwd(codedir)
source('ff_functions_library.R')

#function to compute spatial scales and plot results 
#TODO: add documentation of parameters
calculate_spatial_scales <- function(xs, ys, heading_subsamp = 60, heading_R = 10, speed_dt = 60, fpt_thresh = 600, largest_bin_quantile = 0.999, n_bins = 20, plot = T, short_name = short_name, samprate = samprate, plotdir = '~/Dropbox/cross_species/spatial_scales/'){
  #CALCULATE METRICS 
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  heads <- speeds <- matrix(NA, nrow = n_inds, ncol = n_times)
  
  #indexes to offset for speed computation
  speed_offset <- floor(speed_dt / samprate)
  
  #spatially discretized headings of individuals and speeds of individuals
  for(i in 1:n_inds){
    print(i)
    heads[i,] <- spatial.headings(x=xs[i,],y=ys[i,],R=heading_R,t.idxs=1:n_times,backward = F, fpt.thresh = fpt_thresh, subsample=heading_subsamp)
    idxs_fut <- seq(speed_offset+1,n_times,1)
    idxs_now <- seq(1,n_times-speed_offset,1)
    speeds[i,idxs_now] <- sqrt((xs[i, idxs_fut] - xs[i, idxs_now])^2 + (ys[i, idxs_fut] - ys[i, idxs_now])^2) / (speed_offset * samprate) * 60 #will be in m / min 
  }
  
  #dyadic distances and heading correlations
  dyad_dists <- dyad_dist_changes <- head_corrs <- speed_diffs <- log_speed_diffs <- array(NA, dim = c(n_inds,n_inds,n_times))
  for(i in 1:(n_inds-1)){
    for(j in (i+1):n_inds){
      print(paste(i,j,sep=','))
      dyad_dists[i,j,] <- sqrt((xs[i,]-xs[j,])^2 + (ys[i,]-ys[j,])^2)
      head_corrs[i,j,] <- cos(heads[i,])*cos(heads[j,]) + sin(heads[i,])*sin(heads[j,])
      log_speed_diffs[i,j,] <- abs(log(speeds[i,]) - log(speeds[j,]))
      speed_diffs[i,j,] <- abs(speeds[i,] - speeds[j,])
      idxs_fut <- seq(speed_dt+1,n_times,1)
      idxs_now <- seq(1,n_times-speed_dt,1)
      dyad_dist_changes[i,j,idxs_now] <- dyad_dists[i,j,idxs_fut] - dyad_dists[i,j,idxs_now]
    }
  }
  
  #get distance bins - log spaced from 0 to max distance determined by 99.9% quantile of dyadic distances
  max_dist <- quantile(dyad_dists, largest_bin_quantile, na.rm=T)
  dist_bins <- c(0, 10^seq(0,log(max_dist,10), length.out=n_bins)) #distance bins to use (log spaced)
  
  #hard code distance bins
  dist_bins <- c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192)
  
  #getting means for headings and speed diffs as a function of distance 
  freq_dyad_dists <- heads_mean <- speeds_mean <- change_dyad_dist_mean <- rep(NA, length(dist_bins)-1)
  heads_upper <- speeds_upper <- change_dyad_dist_upper <- rep(NA, length(dist_bins)-1)
  heads_lower <- speeds_lower <- change_dyad_dist_lower <- rep(NA, length(dist_bins)-1)
  for(i in 1:(length(dist_bins)-1)){
    idxs <- which(dyad_dists >= dist_bins[i] & dyad_dists < dist_bins[i+1])
    freq_dyad_dists[i] <- length(idxs)
    
    #mean
    heads_mean[i] <- mean(head_corrs[idxs], na.rm=T)
    speeds_mean[i] <- mean(speed_diffs[idxs], na.rm=T)
    change_dyad_dist_mean[i] <- mean(dyad_dist_changes[idxs], na.rm=T)
    
    #upper quartile
    heads_upper[i] <- quantile(head_corrs[idxs], 0.75, na.rm=T)
    speeds_upper[i] <- quantile(speed_diffs[idxs], 0.75, na.rm=T)
    change_dyad_dist_upper[i] <- quantile(dyad_dist_changes[idxs], 0.75, na.rm=T)
    
    #lower quartile
    heads_lower[i] <- quantile(head_corrs[idxs], 0.25, na.rm=T)
    speeds_lower[i] <- quantile(speed_diffs[idxs], 0.25, na.rm=T)
    change_dyad_dist_lower[i] <- quantile(dyad_dist_changes[idxs], 0.25, na.rm=T)
  }
  
  #PLOT
  
  #middle of the bins for plotting
  mids <- (dist_bins[1:(length(dist_bins)-1)] + dist_bins[2:length(dist_bins)])/2
  
  #SAVE OUTPUT
  out <- list( dist_bins = dist_bins, 
               mids = mids, 
               heads_mean = heads_mean,
               speeds_mean = speeds_mean,
               freq_dyad_dists = freq_dyad_dists,
               change_dyad_dist_mean = change_dyad_dist_mean
  )
  
  #PLOT
  if(plot){
    
    #set up plot
    quartz(width=13, height = 4)
    par(mfrow=c(1,4), mar = c(5,5,1,1))
    non_na_idxs <- which(!is.na(heads_mean))
    
    #don't plot points when there is too little data going into them
    too_little_data_idxs <- which(freq_dyad_dists < 1000)
    heads_mean[too_little_data_idxs] <- NA
    speeds_mean[too_little_data_idxs] <- NA
    change_dyad_dist_mean[too_little_data_idxs] <- NA
    non_na_idxs <- which(!is.na(heads_mean))
  
    #Plot 1: Distribution of log(dyadic distances) between hyenas
    frac_dyad_dists <- freq_dyad_dists / sum(freq_dyad_dists)
    plot(mids,frac_dyad_dists, xlab = 'Distance apart (m)', ylab = 'Probability', pch = 19, col = '#00000066', cex = 1.5, ylim = c(0,max(frac_dyad_dists,na.rm=T)),cex.lab=2, cex.axis=1.5, log='x',type='l',lwd=2, main = short_name)
    abline(v=c(0,1,10,100,1000,10000), lty = 2)
    
    #Plot 2: Mean heading correlation as a function of distance 
    plot(mids,heads_mean, xlab = 'Distance apart (m)', ylab = 'Mean heading correlation', pch = 19, col = '#00000066', cex = 1.5, ylim = c(-1,1),cex.lab=2, cex.axis=1.5, log='x', main = short_name)
    polygon(x=c(mids[non_na_idxs],rev(mids[non_na_idxs])), y = c(heads_upper[non_na_idxs], rev(heads_lower[non_na_idxs])), col = '#00000022', border=NA)
    abline(h=0, lty = 2, lwd = 2)
    abline(v=c(0,1,10,100,1000,10000), lty = 2)
    
    #Plot 3: Mean difference in speed as a function of distance 
    plot(mids,speeds_mean, xlab = 'Distance apart (m)', ylab = 'Mean absolute speed difference (m/min)', pch = 19, col = '#00000066', cex = 2, ylim = c(0,max(speeds_upper,na.rm=T)),cex.lab=2, cex.axis=1.5, log='x', main = short_name)
    polygon(x=c(mids[non_na_idxs],rev(mids[non_na_idxs])), y = c(speeds_upper[non_na_idxs], rev(speeds_lower[non_na_idxs])), col = '#00000022', border=NA)
    abline(h=0, lty = 2, lwd = 2)
    abline(v=c(0,1,10,100,1000,10000), lty = 2)
    
    #Plot 4: Mean change in dyadic distane as a function of original dyadic distance
    plot(mids,change_dyad_dist_mean, xlab = 'Distance apart (m)', ylab = 'Mean speed of divergence (m/min)', pch = 19, col = '#00000066', cex = 2, ylim = c(min(change_dyad_dist_lower,na.rm=T),max(change_dyad_dist_upper,na.rm=T)),cex.lab=2, cex.axis=1.5, log='x', main = short_name)
    polygon(x=c(mids[non_na_idxs],rev(mids[non_na_idxs])), y = c(change_dyad_dist_upper[non_na_idxs], rev(change_dyad_dist_lower[non_na_idxs])), col = '#00000022', border=NA)
    abline(h=0, lty = 2, lwd = 2)
    abline(v=c(0,1,10,100,1000,10000), lty = 2)
    
    dev.copy2pdf(file = paste0(plotdir,'spatial_scales_', short_name,'.pdf'))
  }
  
  invisible(out)
  
}

#MAIN

params <- list(files = files,
               heading_subsamp = heading_subsamp,
               heading_R = heading_R,
               speed_dt = speed_dt,
               fpt_thresh = fpt_thresh,
               largest_bin_quantile = largest_bin_quantile,
               begin_day_remove_sec = begin_day_remove_sec)

#load data and run spatial scales computation, produce plots
scales_data <- list()
for(i in 1:length(files)){

  datafile <- files[i]
  
  print('running analysis on file:')
  print(datafile)
  load(datafile)
  
  #for meerkat data, rename variables to fit with other species
  if(!exists('xs')){
    objects <- ls()
    xidx <- grep('allX', objects, ignore.case=F)
    yidx <- grep('allY', objects, ignore.case=F)
    dayidxidx <- grep('dayIdx', objects, ignore.case=F)
    xname <- objects[xidx]
    yname <- objects[yidx]
    dayidxname <- objects[dayidxidx]
    xs <- get(xname)
    ys <- get(yname)
    dayIdx <- get(dayidxname)
    rm(list=c(xname,yname,dayidxname))
  }
  
  #for coati data, get day break indexes
  if(grepl('coati', short_names[i], ignore.case=T)){
    dayIdx <- c(1,which(diff(ts)>1) + 1)
  }
  
  #remove 10 min from beginning of days, to account for day breaks (and GPS turning on)
  if(exists('dayIdx')){
    for(d in 1:length(dayIdx)){
      if(dayIdx[d] < ncol(xs)){
        xs[, dayIdx[d]:(dayIdx[d]+begin_day_remove_sec)] <- NA
        ys[, dayIdx[d]:(dayIdx[d]+begin_day_remove_sec)] <- NA
      }
    }
  }
  
  samprate <- 1
  if(short_names[i] == 'Hyena_SS2023'){
    samprate <- 30
  }
  
  out <- calculate_spatial_scales(xs = xs, ys = ys, heading_subsamp=heading_subsamp, heading_R = heading_R, short_name = short_names[i], samprate = samprate, n_bins = n_bins, plotdir = plotdir)
  rm(list = c('xs','ys'))
  
  if(exists('dayIdx')){
    rm(list=c('dayIdx'))
  }
  
  #store output
  scales_data[[i]] <- out
  
  
  
}





