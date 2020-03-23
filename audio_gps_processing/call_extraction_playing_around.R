#Playing around with visualizing and detecting calls in wav files

library(seewave)
library(tuneR)
setWavPlayer('/usr/bin/afplay')
library(signal)
library(mFilter)

samp.rate <- 8000
sound.file <- '/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/HM_2017_2/COLLAR2/PROCESSED/HM_HRT_R09_AUDIO_file_6_(2017_08_25-06_44_59)_ASWMUX221110_label.wav'
load('/Users/astrandb/Dropbox/meerkats/data/Kalahari2017/RDATA/calls_all.RData')

true.calls <- calls[which(as.character(calls$dyemark)=='HRT'),]

#min and max times when calls were marked
tmin <- min(true.calls$t0.rec.sec,na.rm=T)
tmax <- max(true.calls$t0.rec.sec,na.rm=T)

#read a short piece of the sound file
idx <-6
t0 <- true.calls$t0.rec.sec[idx]
tf <- true.calls$t0.rec.sec[idx] + 1
wav <- readWave(filename=sound.file,from=t0,to=tf,units='seconds')
#specgram(wav@left,n=128,Fs=8000)
spectro(wav,f=8000,wl=128,ovlp=80)

filt <- remez(15,c(0,200/(samp.rate/2),1000/(samp.rate/2),1),c(0,1,1,0))
filtered <- filter(filt,wav@left)
spectro(filtered,f=8000,wl=128,ovlp=80)
wav2 <- Wave(as.numeric(filtered),samp.rate=samp.rate,bit=16)
asWave(	

w