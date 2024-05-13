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
#(.01 in the code). Store these groups in groups_list, a list of lists whose length
#is equal to n_times.

#Stepping through the groups_list, identify all changes in group membership,
#i.e. consecutive time points when the groups do not match. The algorithm flags
# any change, including instances where individuals disappear or reappear in
# groups due to missing data (but these are later ignored). Store in "changes"
# data frame.
#In a last step, find all changes where the individuals were preserved
#(i.e. individuals did not disappear) and where the number of subgroups either
#increased (a fission) or decreased (a fusion). Store these in a data frame
# called events_detected.
#In cases where the number of subgroups went from 1 to 2 (for a fission) or
#2 to 1 (for a fusion), identify the members of the 2 subgroups as group A and B
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
#     $event_type: "fission" or "fusion"
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
detect_fissions_and_fusions <- function(R_inner, R_outer, xs = xs, ys = ys, ts = ts, names = NULL, verbose = T){
  
  #----Identify subgroups at each point
  if(verbose){print('Identifying subgroups at each point using sticky DBSCAN')}
  #number of inds and times
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  
  #day start indexes
  days <- date(ts)
  day_start_idxs <- c(1, which(diff(days)==1)+1)
  day_start_idxs <- c(day_start_idxs, length(ts)+1)
  
  #need to make sure we don't lose the last day - currently we get an error down the road when we run this though
  #day_start_idxs <- c(day_start_idxs, length(ts)+1)
  
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
    for(t in t_day[1:(length(t_day)-1)]){
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
  
  #get names if not null
  if(!is.null(names)){
    for(i in c(fissions,fusions)){
      changes$group_A[i] <- list(names[changes$group_A_idxs[i][[1]]])
      changes$group_B[i] <- list(names[changes$group_B_idxs[i][[1]]])
    }
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