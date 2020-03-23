#Spatially discretized headings

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
