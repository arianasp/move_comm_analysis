#Detect fissions and fusions using "sticky-DBSCAN" method from Libera et al. 
#Start by defining an adjacency matrix ('together' in the code) of which dyads 
#are "connected" at any moment in time. Dyads are considered to be connected if
#they go within a distance R_inner of one another, and they continue to be 
#connected until they leave a distance R_outer of one another on both ends
#(before and after) of the period where their distance dropped below R_inner. 
#This double threshold makes periods of connectedness more stable by removing 
#the "flicker" that would result from having a single threshold.
#NAs are handled by ending a period of connectedness if an individual has an NA
#at the point immediately before / after (if they were connected after / before 
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
#(labels arbitrary). In cases where there were more subgroups, first ...

#TODO: Check how events are being linked together
#TODO: Something going wrong with NAs - look into it
#TODO: Warnings when can't unambiguoulsy identify groups
#TODO: Generalize day_start_idxs to allow continuous dat

library(lubridate)
library(dbscan)
library(plotrix)

detect_fissions_and_fusions <- function(xs, ys, ts, R_inner, R_outer, names, verbose = T){
  
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
  
  if(verbose){print('Finding fissions and fusions')}
  #go through event times and classify events into types
  changes <- data.frame(tidx = event_times)
  changes$datetime <- ts[changes$tidx]
  changes$event_type <- NA
  changes$n_groups_curr <- unlist(lapply(groups_list, length))[event_times]
  changes$n_groups_next <- unlist(lapply(groups_list, length))[event_times+1]
  changes$inds_curr <- changes$inds_next <- ''
  for(i in 1:nrow(changes)){
    t <- changes$tidx[i]
    groups_curr <- groups_list[[t]]
    groups_next <- groups_list[[t+1]]
    changes$inds_curr[i] <- paste0(which(!is.na(xs[,t])),collapse='_')
    changes$inds_next[i] <- paste0(which(!is.na(xs[,t+1])),collapse='_')
  }
  
  #any times when any individual has missing data cannot be considered as FF events 
  #TODO: relax this assumption!
  changes$event_type[which(changes$n_groups_curr < changes$n_groups_next & changes$inds_curr==changes$inds_next)] <- 'fission'
  changes$event_type[which(changes$n_groups_curr > changes$n_groups_next & changes$inds_curr==changes$inds_next)] <- 'fusion'
  
  #get groups for fissions and fusions
  changes$group_A_idxs <- changes$group_A <- changes$group_B_idxs <- changes$group_B <- list(c(0))
  fissions <- which(changes$event_type=='fission')
  fusions <- which(changes$event_type=='fusion')
  
  #make a column to hold a warning if the subgroups couldn't be unambiguously matched for that event
  changes$warning <- F
  
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
        changes$warning[i] <- T
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
        changes$warning[i] <- T
      }
      
    }
  }
  
  for(i in c(fissions,fusions)){
    changes$group_A[i] <- list(names[changes$group_A_idxs[i][[1]]])
    changes$group_B[i] <- list(names[changes$group_B_idxs[i][[1]]])
  }
  
  #final events dataframe (call it events_detected)
  events_detected <- changes[c(fissions, fusions),c('tidx','datetime','event_type','group_A_idxs','group_B_idxs','group_A','group_B','warning')]
  
  #get number of individuals in each subgroup
  events_detected$n_A <- sapply(events_detected$group_A_idxs, length)
  events_detected$n_B <- sapply(events_detected$group_B_idxs, length)
  
  #sort by time
  events_detected <- events_detected[order(events_detected$tidx),]
  
  if(verbose){print('Done.')}
  
  #return things
  ff_data <- list(events_detected = events_detected, groups_list = groups_list, together = together, changes = changes, R_inner = R_inner, R_outer = R_outer)
  return(ff_data) 
  
}

#Visualize all ff events within a certain time frame
visualize_all_ff_events <- function(ff_data, xs, ys, t_start, t_end, t_step = 5){
  
  #get the time window for plotting events
  t_win <- seq(t_start, t_end)
  
  #get x and y data for all individuals during that time window
  xs_t <- xs[,t_win]
  ys_t <- ys[,t_win]
  
  #get number of individuals
  n_inds <- nrow(xs_t)
  n_times <- ncol(xs_t)
  
  #get together data for the event
  together_t <- ff_data$together[,,t_win]
  
  #get all ff events within that time window
  ff_events <- ff_data$events_detected
  ff_events <- ff_events[which(ff_events$tidx >= t_start & ff_events$tidx <= t_end),]
  inds_all <- unique(c(unlist(ff_events$group_A_idxs), unlist(ff_events$group_B_idxs)))
  if(!sum(inds_all==0)==0){ #TODO: find soure of 0s in inds lists
    inds_all <- inds_all[-which(inds_all==0)]
  }
  print(inds_all)
  inds_all <- inds_all[order(inds_all)]
  
  #colors
  cols <- rep('gray',n_inds)
  
  #get window for plotting
  xmin <- min(xs_t[inds_all,],na.rm=T) - ff_data$R_outer
  xmax <- max(xs_t[inds_all,],na.rm=T) + ff_data$R_outer
  ymin <- min(ys_t[inds_all,],na.rm=T) - ff_data$R_outer
  ymax <- max(ys_t[inds_all,],na.rm=T) + ff_data$R_outer
  
  #make the bottom left = 0,0
  xs_t <- xs_t - xmin
  ys_t <- ys_t - ymin
  
  #create the plot
  quartz()
  for(t in seq(1,n_times,t_step)){
    #set up plot for that time step
    plot(NULL, xlim = c(0, xmax-xmin), ylim = c(0, ymax-ymin), asp = 1,xlab = 'Distance E (m)', ylab = 'Distance N (m)')

    #draw lines connecting individuals that are 'together'
    for(i in 1:(n_inds-1)){
      for(j in (i+1):n_inds){
        if(!is.na(together_t[i,j,t])){
          if(together_t[i,j,t]){
            lines(c(xs_t[i,t], xs_t[j,t]), c(ys_t[i,t], ys_t[j,t]))
          }
        }
      }
    }

    #draw points representing each individual, colored by subgroup
    points(xs_t[,t], ys_t[,t], col = cols, pch = 19)
    
    #plot ff event markers
    #get events to indicate with markers (markers appear for 1 time steps to be visible)
    ff_events_curr <- ff_events[which(abs(ff_events$tidx - (t+t_start)) <= (1*t_step)),]
    
    #if there are any events in the time window, plot each one
    if(nrow(ff_events_curr) > 0){
      for(i in 1:nrow(ff_events_curr)){
        tidx_curr <- ff_events_curr$tidx[i]
        event_type <- ff_events_curr$event_type[i]
        group_A_inds <- ff_events_curr$group_A_idxs[i][[1]]
        group_B_inds <- ff_events_curr$group_B_idxs[i][[1]]
        group_AB_inds <- c(group_A_inds, group_B_inds)
        if(sum(group_AB_inds==0)==0){ #TODO: check why there are still some 0s in the inds lists... should be removed earlier, but this is a hack around it for now
          x_centr <- mean(xs_t[group_AB_inds, tidx_curr-t_start],na.rm=T)
          y_centr <- mean(ys_t[group_AB_inds, tidx_curr-t_start],na.rm=T)
        }
        if(event_type == 'fusion'){
          points(x_centr, y_centr, pch = 1)
        } else{
          points(x_centr, y_centr, pch = 4)
        }
      }
      
    }
    
    
    Sys.sleep(.1)
  }
  
}

visualize_ff_event <- function(ff_data, xs, ys, event_number, t_before = 600, t_after = 600, t_step = 5){
  
  #get the event data from the ff_data object
  event_data <- ff_data$events_detected[event_number,]
  
  #get the time window around the event
  t_win <- seq(event_data$tidx-t_before, event_data$tidx + t_after)
  
  #get x and y data for all individuals during that time window
  xs_t <- xs[,t_win]
  ys_t <- ys[,t_win]
  
  #get number of individuals
  n_inds <- nrow(xs_t)
  n_times <- ncol(xs_t)
  
  #get together data for the event
  together_t <- ff_data$together[,,t_win]
  
  #get event type (fission on fusion)
  event_type <- event_data$event_type
  
  #inds in group A and group B
  inds_A <- event_data$group_A_idxs[[1]]
  inds_B <- event_data$group_B_idxs[[1]]
  inds_all <- c(inds_A, inds_B)
  
  #set colors based on group
  cols <- rep('gray', n_inds)
  cols_inner <- rep('#00000005',n_inds)
  cols_outer <- rep('#00000003',n_inds)
  cols[inds_A] <- 'red'
  cols[inds_B] <- 'blue'
  cols_inner[inds_A] <- '#FF000005'
  cols_inner[inds_B] <- '#0000FF05'
  cols_outer[inds_A] <- '#FF000003'
  cols_outer[inds_B] <- '#0000FF03'
  
  
  #get window for plotting
  xmin <- min(xs_t[inds_all,],na.rm=T) - ff_data$R_outer
  xmax <- max(xs_t[inds_all,],na.rm=T) + ff_data$R_outer
  ymin <- min(ys_t[inds_all,],na.rm=T) - ff_data$R_outer
  ymax <- max(ys_t[inds_all,],na.rm=T) + ff_data$R_outer
  
  #make the bottom left = 0,0
  xs_t <- xs_t - xmin
  ys_t <- ys_t - ymin
  
  #create the plot
  quartz()
  for(t in seq(1,n_times,t_step)){
    #set up plot for that time step
    plot(NULL, xlim = c(0, xmax-xmin), ylim = c(0, ymax-ymin), asp = 1,xlab = 'Distance E (m)', ylab = 'Distance N (m)', main = paste(event_data$event_type))
    
    #draw lines connecting individuals that are 'together'
    for(i in 1:(n_inds-1)){
      for(j in (i+1):n_inds){
        if(!is.na(together_t[i,j,t])){
          if(together_t[i,j,t]){
            lines(c(xs_t[i,t], xs_t[j,t]), c(ys_t[i,t], ys_t[j,t]))
          }
        }
      }
    }
    
    #draw points representing each individual, colored by subgroup
    points(xs_t[,t], ys_t[,t], col = cols, pch = 19)
    Sys.sleep(.1)
  }
  
}


