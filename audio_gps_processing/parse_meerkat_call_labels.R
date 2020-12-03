library(docstring)

hybrid_to_standard_format <- function(label){
  #' Parses a single hybrid call to the standard format for hybrids
  #' hyb a+b where a and b are the component call types in alphabetical order
  #' @param label: a hybrid call label (string)
  #' @return parsed.label: the call label in standard form
  
  name.split <- unlist(strsplit(label,' '))
  hybseq <- name.split[which(grepl('\\+', name.split))]
  splitseq <- unlist(strsplit(hybseq, '\\+'))
  sorted <- sort(splitseq)
  collapsed <- paste(sorted,collapse='+')
  parsed.label <- paste('hyb', collapsed, sep = ' ')
  
  return(parsed.label)
  
}

parse_meerkat_call_labels <- function(call.labels){
  #' Parses call lables into a specified set of categories
  #' 
  #' For callType columnm, hybrids are retained.
  #' For simpleType, hybrids are removed, following the rules below, resulting in the following 'simple' categories:
  #' cc, cchyb (hybrid cc), agg, soc, chat, mo, ld, chat, al, oth (other)
  #' Rules for 'simple' method:
  #' a "hierarchy" of call priority (when hybrids exist, they get classified with the higher up category)
  #' hierarchy: cc, al, ld, mo, chat, soc, agg, s, lc, oth
  
  #' @param call.labels: vector of character strings giving call labels (should be input from file 2017_ALL_CALLS_SYNCHED.csv or 2019_ALL_CALLS_SYNCHED.csv)

  #' @return parsed.calls: data frame containing columns: 
  #' $callType or callSimple (parsed call type, depending on if standard or simple method selected)
  #' $isCall (1 if it is a call, 0 if not (e.g. synch), NA if unsure ('#' indicated))
  #' $focal (1 if focal, 0 if nonfocal, -1 if unsure, NA if not a call or unsure call)
  #' $noisy (1 if noisy, 0 if not noisy, NA if not a call or unsure call)
  #' $unsureType (1 if unsure of the type ('?' indicated), 0 if not, NA if not a call or unsure call)
  #' hybrid (1 if the call is a hybrid, 0 if not, NA if not a call)
  
  #convert labels to lower case
  call.labels <- tolower(as.character(call.labels))
  
  #remove anything in parentheses
  call.labels <- gsub("\\([^\\]]*\\)", "", call.labels, perl=TRUE)
  
  #get number of calls
  n <- length(call.labels)
  
  #create data frame to hold results. First assume all calls are other
  calls.parsed <- data.frame(callType = rep('oth', n), 
                             callSimple = rep('oth', n), 
                             isCall = rep(NA, n),
                             focal = rep(NA, n),
                             noisy = rep(NA, n),
                             unsureType = rep(NA, n), 
                             hybrid = rep(NA,n), 
                             stringsAsFactors = F)
  
  #deal with the non hybrid calls
  call.types.simple <- c('s','lc','agg','soc','chat','mo','ld','al','cc','ukn','oth')
  calls.parsed$callSimple[which(grepl('s', call.labels,ignore.case=T))] <- 's'
  calls.parsed$callSimple[which(grepl('lc',call.labels,ignore.case=T))] <- 'lc'
  calls.parsed$callSimple[which(grepl('lost',call.labels,ignore.case=T))] <- 'lc'
  calls.parsed$callSimple[which(grepl('ag',call.labels,ignore.case=T))] <- 'agg'
  calls.parsed$callSimple[which(grepl('so',call.labels,ignore.case=T))] <- 'soc'
  calls.parsed$callSimple[which(grepl('chat',call.labels,ignore.case=T))] <- 'chat'
  calls.parsed$callSimple[which(grepl('mo',call.labels,ignore.case=T))] <- 'mo'
  calls.parsed$callSimple[which(grepl('ld',call.labels,ignore.case=T))] <- 'ld'
  calls.parsed$callSimple[which(grepl('lead',call.labels,ignore.case=T))] <- 'ld'
  calls.parsed$callSimple[which(grepl('al',call.labels,ignore.case=T))] <- 'al'
  calls.parsed$callSimple[which(grepl('cc',call.labels,ignore.case=T))] <- 'cc'
  calls.parsed$callSimple[which(call.labels == 'c')] <- 'cc'
  calls.parsed$callSimple[which(grepl('Marker',call.labels,ignore.case=T))] <- 'cc'
  calls.parsed$callSimple[which(grepl('uk',call.labels,ignore.case=T))] <- 'ukn'
  calls.parsed$callSimple[which(grepl('unk',call.labels,ignore.case=T))] <- 'ukn'
  
  #deal with synch calls, beeps, skips, and various other non-call things
  calls.parsed$callSimple[which(grepl('syn',call.labels,ignore.case=T))] <- 'synch'
  calls.parsed$callSimple[which(grepl('syc',call.labels,ignore.case=T))] <- 'synch'
  calls.parsed$callSimple[which(grepl('bee',call.labels,ignore.case=T))] <- 'beep'
  calls.parsed$callSimple[which(grepl('skipon',call.labels,ignore.case=T))] <- 'skipon'
  calls.parsed$callSimple[which(grepl('skip on',call.labels,ignore.case=T))] <- 'skipon'
  calls.parsed$callSimple[which(grepl('skipoff',call.labels,ignore.case=T))] <- 'skipoff'
  calls.parsed$callSimple[which(grepl('skip off',call.labels,ignore.case=T))] <- 'skipoff'
  calls.parsed$callSimple[which(grepl('start',call.labels,ignore.case=T))] <- 'start'
  calls.parsed$callSimple[which(grepl('stop',call.labels,ignore.case=T))] <- 'stop'
  calls.parsed$callSimple[which(grepl('end',call.labels,ignore.case=T))] <- 'stop'
  calls.parsed$callSimple[which(grepl('dig',call.labels,ignore.case=T))] <- 'dig'
  calls.parsed$callSimple[which(grepl('eat',call.labels,ignore.case=T))] <- 'eat'
  calls.parsed$callSimple[which(grepl('chew',call.labels,ignore.case=T))] <- 'eat'
  calls.parsed$callSimple[which(grepl('pause',call.labels,ignore.case=T))] <- 'pause'
  calls.parsed$callSimple[which(grepl('oor',call.labels,ignore.case=T))] <- 'oor'
  calls.parsed$callSimple[which(grepl('bir',call.labels,ignore.case=T))] <- 'bir'
  calls.parsed$callSimple[which(grepl('outrange',call.labels,ignore.case=T))] <- 'oor'
  calls.parsed$callSimple[which(grepl('inrange',call.labels,ignore.case=T))] <- 'bir'
  
  #create column for the simple call types
  calls.parsed$callType <- calls.parsed$callSimple
  
  #then find the hybrids and deal with them
  #all hybrids (including sq and fu) receive a label hyb a+b where a and b are component calls in alphabetical order
  hyb.idxs <- unique(c(which(grepl('hyb', call.labels)), which(grepl('fu', call.labels)), which(grepl('sq', call.labels))))
  hyb.labels <- sapply(call.labels[hyb.idxs], FUN = hybrid_to_standard_format)
  calls.parsed$callType[hyb.idxs] <- hyb.labels
  
  #get whether it's a call
  calls.parsed$isCall <- 0
  calls.parsed$isCall[which(calls.parsed$callSimple %in% call.types.simple)] <- 1
  calls.parsed$isCall[which(grepl('\\#', call.labels, ignore.case=T))] <- NA
  
  #get focal/nonfocal info
  calls.parsed$focal <- 1
  calls.parsed$focal[which(grepl('nf', call.labels, ignore.case=T))] <- 0
  calls.parsed$focal[which(grepl('non', call.labels, ignore.case=T))] <- 0
  calls.parsed$focal[which(grepl('\\*', call.labels, ignore.case=T))] <- -1
  calls.parsed$focal[which(calls.parsed$isCall==0 | is.na(calls.parsed$isCall))] <- NA
  
  #get noisy info
  calls.parsed$noisy <- 0
  calls.parsed$noisy[which(grepl('x', call.labels, ignore.case=T))] <- 1
  calls.parsed$noisy[which(calls.parsed$isCall==0 | is.na(calls.parsed$isCall))] <- NA
  
  #get whether unsure about type
  calls.parsed$unsureType <- 0
  calls.parsed$unsureType[which(grepl('\\?', call.labels, ignore.case=T))] <- 1
  calls.parsed$unsureType[which(calls.parsed$isCall==0 | is.na(calls.parsed$isCall))] <- NA

  #get whether it's a hybrid
  calls.parsed$hybrid <- 0
  calls.parsed$hybrid[hyb.idxs] <- 1
  calls.parsed$hybrid[which(calls.parsed$isCall==0 | is.na(calls.parsed$isCall))] <- NA
  
  return(calls.parsed)
  
}



