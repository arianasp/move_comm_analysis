#This script identifies possible pairs of calls labeled as 'focal' on different devices that might 'match' (i.e. actually be the same call)
#It then takes you through these matches and you can manually identify whether they are, in fact, the same
#After some setup, the script will begin taking you through pairs of calls and you will be able to interactively label them as matches or not (or unsure)
#By default, the script just plays you the calls in sequence, however you can also choose to generate a spectrogram for further investigation
#You can quit at any time and your labels will be saved
#At the moment I haven't implemented any way to go back and correct things if you make a mistake, so just try not to make mistakes!
#But in case you notice you've made a mistake, note down the two id numbers shown when that pair of calls is played and save them somewhere - this will make it easy to go back and fix it later

#---------------- INSTRUCTIONS FOR HOW TO RUN ------------------
#In the next section, enter in the appropriate parameters for each line marked "INPUT"
#After you have completed this, save the script
#To run, first type into the R console: setwd(DIR) where DIR is the directory in which you've stored the code
#Then, type into the R consolte: source('identify_possible_misidentified_focal_calls.R') to run the script
#Note that the script can't be nicely run by shift-enter because of the interactive bits, so that's why you have to use the 'source' way instead.

#----------------PARAMETERS--------------------

#TO RUN: Input the following fields

#INPUT: your name
name <- 'ari'

#INPUT: which year's data to use
year <- 2019

#INPUT: which meerkat group
group <- 'HM'

#INPUT: directory where call match data are/will be stored
match.save.dir <- '~/Dropbox/meerkats/meerkats_shared/ari'

#INPUT: directories where data is stored
label.dir <- '~/Dropbox/meerkats/meerkats_shared/data'
audio.dir <- paste('/Volumes/Public/MEERKAT_RAW_DATA/', year, sep='')
gps.dir <- paste('~/Dropbox/meerkats/meerkats_shared/data')

#INPUT: wav player - this might be mac/pc difference but /usr/bin/afplay works for me
wav.player <- '/usr/bin/afplay'

#LEAVE THE PARAMETERS BELOW AS IS

#time threshold for how close together two calls have to be to be considered possibly the same
time.thresh <- 0.3

#time threshold for how similar in duration two calls have to be to be flagged as possibly the same
dur.thresh <- 0.1

#padding to play on either side of a call
pad <- 0.05

#maximum distance apart to check (to be conservative, set this to 20m)
max.dist.apart <- 20

#---------------SETUP---------------------

#source functions
source('play_compare_calls.R')

#label file name
label.file <- paste(label.dir,'/',year, '_ALL_CALLS_SYNCHED.csv', sep = '')

#gps file name
gps.file <- paste(gps.dir,'/',group,'_COORDINATES_',year,'_sessions.RData',sep='')

#individual info file name
ind.info.file <- paste(gps.dir , '/', group, '_INDIVIDUAL_INFO_', year, '.txt', sep = '')

#whether to load matches or generate them
load.matches <- T
loadfromfile <- readline('Would you like to load matches table from a previous run of the script (y) or start over (n)? ')
if(loadfromfile=='n'){
  print('Warning, any manual checks of matches that you had saved previously will be erased if you overwrite a new file - only do this if this is the first time youre running this script!')
  ok <- readline('Is that all OK? (y/n) ')
  if(ok == 'y'){
    load.matches <- F
  }
}

print('setting up, please wait...')

#set wave player
setWavPlayer(wav.player)

#set working directory to audio directory
setwd(audio.dir)

#set file where matches are or will be saved, based on name of labeler and year and directory
match.save.file <- paste(match.save.dir, '/', 'matches_', year, '_', name, '.RData', sep = '')

#read in labels file if needed
if(!load.matches){
  calls <- read.csv(label.file,sep='\t', stringsAsFactors=F)
}

#get paths to all audio files in that year
if(!load.matches){
  all.audio.files <- list.files(path = audio.dir, recursive = T, pattern = '*.wav', ignore.case = T)
}

#convert t0 and tf to number
if(!load.matches){
  calls$t0num <- as.numeric(as.POSIXlt(calls$t0GPS))

  #give each call a unique identifier
  calls$unique.id <- seq(1, nrow(calls), 1)
}

#get gps time index for


#if needed, create matches data frame

if(!load.matches){
  #loop over filenames and flag possible matches
  matches <- data.frame()

  files.all <- unique(calls$fileName)

  for(i in 1:length(files.all)){

    #separate out one file
    calls.file <- calls[which(calls$fileName == files.all[i] & calls$duration > 0),]
    calls.nonfile <- calls[which(calls$fileName != files.all[i] & calls$duration > 0),]

    #for each row in that file, search for matches
    for(j in 1:nrow(calls.file)){

      t0.call <- calls.file$t0num[j]
      dur.call <- calls.file$duration[j]

      #only look at calls with duration > 0 sec (no 'start' markers included)
      if(dur.call > 0){

        #calculate time diff between calls from other files and this call
        dt <- abs(calls.nonfile$t0num - t0.call)

        #calculate difference in duration between calls from other files and this call
        ddur <- abs(calls.nonfile$duration - dur.call)

        #find potential matches
        idx.match <- which(dt <= time.thresh & ddur <= dur.thresh)

        #get number of matches
        n.matches <- length(idx.match)

        #if some matches found, add information to data frame of possible matches
        if(n.matches>0){

          match.rows <- data.frame(unique.id.a = rep(calls.file$unique.id[j], n.matches), unique.id.b = calls.nonfile$unique.id[idx.match])
          matches <- rbind(matches, match.rows)

        }
      }


    }

  }

  #first relabel ids of the two events by so that a is min and b is max (arbitrary, but helps with removing duplicates)
  matches$min.unique.id <- apply(cbind(matches$unique.id.a, matches$unique.id.b),1,min)
  matches$max.unique.id <- apply(cbind(matches$unique.id.a, matches$unique.id.b),1,max)
  matches$unique.id.a <- matches$min.unique.id
  matches$unique.id.b <- matches$max.unique.id
  matches <- matches[,c('unique.id.a','unique.id.b')]

  #remove duplicates (since all overlap events have gotten picked up twice, once on each file)
  dups <- duplicated(cbind(matches$unique.id.a,matches$unique.id.b))
  matches <- matches[!dups,]

  #add some other useful columns
  matches$file.a <- calls$fileName[match(matches$unique.id.a, calls$unique.id)]
  matches$file.b <- calls$fileName[match(matches$unique.id.b, calls$unique.id)]
  matches$t0.rec.a <- calls$t0File[match(matches$unique.id.a, calls$unique.id)]
  matches$t0.rec.b <- calls$t0File[match(matches$unique.id.b, calls$unique.id)]
  matches$dur.a <- calls$duration[match(matches$unique.id.a, calls$unique.id)]
  matches$dur.b <- calls$duration[match(matches$unique.id.b, calls$unique.id)]
  matches$dt <- abs(calls$t0num[match(matches$unique.id.a, calls$unique.id)] - calls$t0num[match(matches$unique.id.b, calls$unique.id)])
  matches$type.a <- calls$callType[match(matches$unique.id.a, calls$unique.id)]
  matches$type.b <- calls$callType[match(matches$unique.id.b, calls$unique.id)]
  matches$nonfoc.a <- calls$nonFocal[match(matches$unique.id.a, calls$unique.id)]
  matches$nonfoc.b <- calls$nonFocal[match(matches$unique.id.b, calls$unique.id)]
  matches$unsurefoc.a <- calls$unsureFocal[match(matches$unique.id.a, calls$unique.id)]
  matches$unsurefoc.b <- calls$unsureFocal[match(matches$unique.id.b, calls$unique.id)]
  matches$ind.a <- calls$ind[match(matches$unique.id.a, calls$unique.id)]
  matches$ind.b <- calls$ind[match(matches$unique.id.b, calls$unique.id)]


  #get rid of synchs and beeps
  matches <- matches[which(!(matches$type.a == 'synch')),]
  matches <- matches[which(!(matches$type.b == 'synch')),]
}

#load gps data and add distance info to the matches table
if(!load.matches){

  load(gps.file)
  ind.info <- read.table(ind.info.file, header=T, sep = '\t')

  #add columns to matches table with index in gps data matrices
  matches$ind.idx.a <- match(matches$ind.a,ind.info$code)
  matches$ind.idx.b <- match(matches$ind.b,ind.info$code)

  #get time index
  matches$tGPS <- calls$t0GPS_UTC[match(matches$unique.id.a, calls$unique.id)]
  matches$date <- calls$date[match(matches$unique.id.a, calls$unique.id)]
  matches$date.idx <- match(matches$date, dates)
  matches$t.idx <- calls$t0_idx_dayTimeline[match(matches$unique.id.a, calls$unique.id)] + dayIdx[matches$date.idx] -1

  #get x and y coordinates
  matches$x.a <- allX[cbind(matches$ind.idx.a, matches$t.idx)]
  matches$y.a <- allY[cbind(matches$ind.idx.a, matches$t.idx)]
  matches$x.b <- allX[cbind(matches$ind.idx.b, matches$t.idx)]
  matches$y.b <- allY[cbind(matches$ind.idx.b, matches$t.idx)]

  #get distance between a and b
  matches$dist.apart <- sqrt( (matches$x.a - matches$x.b)^2 + (matches$y.a - matches$y.b)^2 )

}

#make a table that matches label files to wav files for easy access
if(!load.matches){
  labels.to.wav.files <- data.frame(label.file = unique(c(matches$file.a, matches$file.b)))
  labels.to.wav.files$label.file <- as.character(labels.to.wav.files$label.file)
  labels.to.wav.files$wav.file <- NA

  for(idx in 1:nrow(labels.to.wav.files)){
    label.file <- labels.to.wav.files$label.file[idx]
    wav.file <- NULL
    found <- F
    for(i in 1:length(all.audio.files)){

      ismatch <- grepl(escapeRegex(gsub('.wav','',basename(all.audio.files[i]),ignore.case=T)), label.file)
      if(ismatch){
        wav.file <- all.audio.files[i]
        found <- T
        labels.to.wav.files$wav.file[idx] <- wav.file
      }
    }
    if(!found){
      warning(paste('matching file not found for label file',label.file))
    }
  }
}

#add a column to indicate distance between meerkats

#find indices of matches to look at (where both calls are labeled as focal and distance < max.dist.apart or NA)
if(!load.matches){
  idxs <- which(matches$nonfoc.a == 0
                & matches$nonfoc.b == 0
                & matches$unsurefoc.a == 0
                & matches$unsurefoc.b == 0
                & (is.na(matches$dist.apart) | (matches$dist.apart < 20)))
}

#if running for first time, initialize manual label column
if(!load.matches){
  matches$manual.check.result <- NA
}

if(load.matches){
  load(match.save.file)
}

print(paste('done setting up. there are ', length(idxs), ' matched calls to look at and ', sum(!is.na(matches$manual.check.result)), ' have been already manually checked.', sep = ''))

#--------------------MAIN--------------------------

print('ready to label some calls? (y/n)')
ready <- readline()

if(ready == 'y'){
  #Loop over indexes we are checking to listen / look at spectrograms and give a manual label
  for(i in 1:length(idxs)){

    row <- idxs[i]

    #if the call pair was already checked, skip it
    if(!is.na(matches$manual.check.result[row])){
      next
    } else{

      print(paste('now viewing calls: ', matches$unique.id.a[row], ' AND ' ,matches$unique.id.b[row], sep = ''))


      #compare calls and store output in matches data frame
      matches$manual.check.result[row] <- play.compare.calls(matches, labels.to.wav.files, row)
      if(is.na(matches$manual.check.result[row])){
        break
      }
    }
  }

  #save output

  cat('Would you like to save the labels you created? (y/n)')
  saveorno <- readline()

  if(saveorno == 'y'){
    save(file=match.save.file,list=c('matches', 'calls', 'time.thresh','dur.thresh','pad','year','labels.to.wav.files','idxs','all.audio.files'))
    print('successfully saved to: ')
    print(match.save.file)
  } else{
    print('WARNING: table was not saved')
  }
}

