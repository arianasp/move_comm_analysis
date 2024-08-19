#Playing with analysis of call sequences

library(lubridate)
library(fields)
library(viridisLite)
library(ks)

#----FUNCS---
#get transition matrix
#transition_matrix[i,j] is the number of times call type i transitions to call type j
get_transition_matrix <- function(preds){
  call_types <- names(sort(table(preds$label), decreasing = T))
  n_calls <- length(call_types)
  transition_matrix <- matrix(NA, nrow = n_calls, ncol = n_calls)
  for(i in 1:n_calls){
    for(j in 1:n_calls){
      idxs <- which(preds$label == call_types[i] & preds$next_label == call_types[j])
      transition_matrix[i,j] <- length(idxs)
    }
  }
  out <- list()
  out$call_types <- call_types
  out$transition_matrix <- transition_matrix
  return(out)
}

#randomize call order by shuffling labels
#also update next_call to reflect the new updated order
#return a data frame preds_rand
randomize_call_order <- function(preds){
  filenames <- unique(preds$filename)
  preds_rand <- preds
  preds_rand$label <- NA
  for(f in 1:length(filenames)){
    curr_idxs <- which(preds$filename == filenames[f])
    if(length(curr_idxs)>1){
      preds_rand$label[curr_idxs] <- sample(preds$label[curr_idxs])
      next_idxs <- curr_idxs[2:length(curr_idxs)]
      preds_rand$next_label[curr_idxs] <- c(preds_rand$label[next_idxs], NA)
    }
  }
  return(preds_rand)
}


#directory holding label files
indir <- '~/EAS_shared/meerkat/working/processed/acoustic/animal2vec_predictions/large_model/2017/csv'

#----MAIN---

#GET DATA
setwd(indir)

files <- list.files(pattern = 'csv')
files_primary <- files[-grep('_secondary_predictions', files)]

#load prediction files and save to a big table
preds_all <- data.frame()
for(i in 1:length(files_primary)){

  #get current file and get name without extension
  file_curr <- files_primary[i]
  file_no_ext <- gsub(pattern = "\\.csv$", "", file_curr)
  
  #get group and id from filename
  group <- strsplit(file_no_ext,'_')[[1]][1]
  id <- strsplit(file_no_ext,'_')[[1]][2]
  
  #get predictions from the current file
  preds_curr <- read.csv(file = file_curr, header=T, sep = '\t')
  
  #create a data table with individual id, filename, start time in sec, end time in sec
  preds_curr_formatted <- data.frame(filename = rep(file_no_ext, nrow(preds_curr)),
                                     group = rep(group, nrow(preds_curr)),
                                     id = rep(id, nrow(preds_curr)),
                                     label = preds_curr$Name,
                                     start_time = as.numeric(hms(preds_curr$Start)),
                                     duration = as.numeric(hms(preds_curr$Duration)))
  preds_curr_formatted$end_time <- preds_curr_formatted$start_time + preds_curr_formatted$duration
  
  #sort by start_time (though they should already be sorted)
  preds_curr_formatted <- preds_curr_formatted[order(preds_curr_formatted$start_time),]
  
  preds_all <- rbind(preds_all, preds_curr_formatted)
}

#remove non focal calls and non calls
nonfoc_rows <- grep('nf',preds_all$label)
noncall_rows <- which(preds_all$label %in% c('beep','synch','eating'))
rows_to_remove <- unique(c(nonfoc_rows, noncall_rows))
if(length(rows_to_remove)>0){
  preds <- preds_all[-rows_to_remove,]
} 

#get overall frequency of call types
sort(table(preds$label), decreasing = T)

#get next call and inter-call interval
filenames <- unique(preds$filename)
preds$next_label <- NA
preds$inter_call_interval <- NA
for(f in 1:length(filenames)){
  curr_idxs <- which(preds$filename == filenames[f])
  if(length(curr_idxs)>1){
    next_idxs <- curr_idxs[2:length(curr_idxs)]
    preds$next_label[curr_idxs] <- c(preds$label[next_idxs], NA)
    preds$inter_call_interval[curr_idxs[1:(length(curr_idxs)-1)]] <- c(preds$start_time[next_idxs] - preds$start_time[curr_idxs[1:(length(curr_idxs)-1)]])
  }
}



out <- get_transition_matrix(preds)
call_types <- out$call_types
transition_matrix <- out$transition_matrix
n_calls <- length(call_types)

#noramlize transition matrix
#transition_matrix_norm[i,j] is hte probability that give you've just made call i you transition to call j
transition_matrix_norm <- transition_matrix
for(j in 1:n_calls){
  transition_matrix_norm[,j] <- transition_matrix[,j] / sum(transition_matrix[,j], na.rm=T)
}

#get randomized call transition matrices
n_rands <- 100
transition_matrices_rand <- array(NA, dim = c(n_calls, n_calls, n_rands))
for(i in 1:n_rands){
  print(i)
  preds_rand <- randomize_call_order(preds)
  transition_matrices_rand[,,i] <- get_transition_matrix(preds_rand)$transition_matrix
}


#----PLOTS---
cols <- viridis(n_calls)
pdf(file = '~/call_seqs_analyses.pdf',width = 14, height = 5)
par(mfrow=c(1,3))

#CALL TYPE FREQUENCY
freqs <- sort(table(preds$label), decreasing = T)
fracs <- freqs / sum(freqs) * 100
barplot(fracs, names.arg = names(freqs), col = cols, border = cols, xlab = 'Call type', ylab = 'Percentage')

#transition matrix (normalized)
#image.plot(t(transition_matrix_norm), col = viridis(256), xaxt = 'n', yaxt = 'n', xlab = 'First call', ylab = 'P(second call | first call)')
#axis(1, at = seq(0,1,length.out=n_calls), labels = call_types)
#axis(2, at = seq(0,1,length.out=n_calls), labels = call_types)

#----under and over representation of call transitions
transition_matrix_rand_mean <- apply(transition_matrices_rand, MARGIN = c(1,2), FUN = mean)
transition_matrix_rand_upper <- apply(transition_matrices_rand, MARGIN = c(1,2), FUN = function(x){return(quantile(x, 0.975, na.rm=T))})
transition_matrix_rand_lower <- apply(transition_matrices_rand, MARGIN = c(1,2), FUN = function(x){return(quantile(x, 0.025, na.rm=T))})
log_ratios <- log(transition_matrix / transition_matrix_rand_mean)
colpal <- colorRampPalette(c('purple','white','red'))
image.plot(t(log_ratios), col = colpal(256), xaxt = 'n', yaxt = 'n', xlab = 'First call', ylab = 'Second call',zlim=c(-2.5,2.5))
axis(1, at = seq(0,1,length.out=n_calls), labels = call_types)
axis(2, at = seq(0,1,length.out=n_calls), labels = call_types)

#INTER-CALL INTERVALS
xmax <- 100
ymax <- 1.2
types_to_plot <- call_types
cols <- viridis(length(types_to_plot))
ltys <- rep(seq(1,3),3)
plot(NULL, xlim = c(.1, 5), ylim = c(0,ymax), xlab = 'Inter-call interval (sec)', ylab = 'Density')
for(i in 1:length(types_to_plot)){
  call_type_for_interval <- types_to_plot[i]
  intervals <- preds$inter_call_interval[which(preds$label == call_type_for_interval & preds$next_label == call_type_for_interval)]
  kdes <- kde(intervals, h = .2, gridsize = 10000, xmin = 0, xmax = xmax)
  lines(kdes$eval.points, kdes$estimate, col = cols[i], lwd = 3, lty = ltys[i])
}
legend('topright',legend = types_to_plot, col = cols, lwd = 3, lty = ltys)

dev.off()


#visualize a time series of calls
pdf(file = '~/call_timeseries_example.pdf', width = 14, height = 4)
file <- filenames[2]
preds_curr <- preds[which(preds$filename == file),]

plot(NULL, xlim=c(0,max(preds_curr$end_time, na.rm=T))/60, ylim = c(0,n_calls), xlab = 'Time (min)' , ylab = 'Call type', yaxt = 'n')
abline(h = 1:n_calls, col = '#00000022')
for(i in 1:n_calls){
  preds_of_type <- preds_curr[which(preds_curr$label == call_types[i]),]
  arrows(preds_of_type$start_time/60, rep(i-1, nrow(preds_of_type)), preds_of_type$start_time/60, rep(i, nrow(preds_of_type)) , col = cols[i], length = 0)
}
axis(side = 2, at = seq(0.5,n_calls,1), labels = call_types, las = 2)
dev.off()


