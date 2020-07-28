#--------------FUNCTIONS FOR FOCAL / NON-FOCAL IDENTIFICATION ----------

#--------------LIBRARIES------------------
library(Hmisc)
library(tuneR)
library(seewave)
library(lubridate)
library(signal)
library(oce)

#--------------FUNCTIONS------------------------
#extract and line up spectrograms of matched calls for visual and audio comparison
play.compare.calls <- function(matches, labels.to.wav.files, i, pad=.1){

  play.again <- T
  while(play.again){
    file.a <- labels.to.wav.files$wav.file[which(labels.to.wav.files$label.file==matches$file.a[i])]
    file.b <- labels.to.wav.files$wav.file[which(labels.to.wav.files$label.file==matches$file.b[i])]
    t0.a <- as.numeric(seconds(hms(matches$t0.rec.a[i]))) - pad
    tf.a <- t0.a + matches$dur.a[i] + 2*pad
    t0.b <- as.numeric(seconds(hms(matches$t0.rec.b[i]))) - pad
    tf.b <- t0.b + matches$dur.b[i] + 2*pad

    #make spectrograms the same size for easier comparison
    tot.len <- max(tf.a - t0.a, tf.b - t0.b)

    wav.a <- readWave(file.a, from = t0.a, to = t0.a + tot.len, units='seconds')
    wav.b <- readWave(file.b, from = t0.b, to = t0.b + tot.len, units='seconds')
    specwin <- 256

    if(wav.a@samp.rate > 8000){
      wav.a <- downsample(wav.a, 8000)
    }
    if(wav.b@samp.rate > 8000){
      wav.b <- downsample(wav.b, 8000)
    }

    if(length(dev.list() > 20)){
      dev.off(dev.list()["RStudioGD"])
    }

    par(mfrow=c(1,2), mar = c(1,1,1,1))

    #first spectrogram
    spec.a <- spectro(wav.a, wl = specwin, ovlp=95, flog=T, plot = F, wn = 'blackman')
    image.plot(t(spec.a$amp), x = spec.a$time, y = spec.a$freq, horizontal = T, xlab = 'n', ylab = 'n', col = viridis(1024))
    start.frac.a <- pad / tot.len
    end.frac.a <- (tf.a - t0.a - pad) / tot.len
    start.spec.a <- start.frac.a * max(spec.a$time)
    end.spec.a <- end.frac.a * max(spec.a$time)
    abline(v = c(start.spec.a, end.spec.a), lwd = 2, col = 'black', lty = 2)

    #second spectrogram
    spec.b <- spectro(wav.b, wl = specwin, ovlp=95, flog=T, plot = F, wn = 'blackman')
    image.plot(t(spec.b$amp), x = spec.b$time, y = spec.b$freq, horizontal = T, xlab = 'n', ylab = 'n', col = viridis(1024))
    start.frac.b <- pad / tot.len
    end.frac.b <- (tf.b - t0.b - pad) / tot.len
    start.spec.b <- start.frac.b * max(spec.b$time)
    end.spec.b <- end.frac.b * max(spec.b$time)
    abline(v = c(start.spec.b, end.spec.b), lwd = 2, col = 'black', lty = 2)

    play(wav.a)
    play(wav.b)

    cat('Do they look / sound the same? (y = yes, n = no, u = unsure, a = play again, s = draw spectrogram, q = save and quit)')
    user.input <- readline()

    if(user.input %in% c('y','n','u')){
      play.again <- F

      if(user.input == 'y'){
        cat('Which one do you think is the focal? (1 = left/first, 2 = right/second, u = unknown) You can press a to replay ')
        whichfoc <- readline()
        while(whichfoc == 'a'){
          play(wav.a)
          play(wav.b)
          whichfoc <- readline()
        }
        out <- paste(user.input, whichfoc, sep = '_')
        return(out)
      } else{


        return(user.input)
      }
    } else if(user.input == 'q'){
      return(NA)
    }

  }

}

