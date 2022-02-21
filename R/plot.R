plot.Noblesse <- function(x,run=1,show.zeros=TRUE){
    method <- mdf(dat[[run]]$mdf$name)
    signal <- x[[run]]$signal
    nz <- ifelse(show.zeros,2,0)
    nm <- length(method$codes)
    nc <- ceiling(sqrt(nm+nz))
    nr <- ceiling((nm+nz)/nc)
    oldpar <- par(mfrow=c(nr,nc),mar=c(4,3.5,1,1),oma=c(0,0,3,0),mgp=c(2,1,0))
    on.exit(par(oldpar))
    plot(x=signal[['1']][,3],y=signal[['1']][,1],
         xlab='t(s)',ylab='Far',pch=16,main='1')
    plot(x=signal[['1']][,3],y=signal[['1']][,2],
         xlab='t(s)',ylab='IC0',pch=16,main='1')
    for (i in 1:nm){
        code <- method$codes[i]
        column <- method$columns[i]
        label <- method$labels[i]
        plot(x=signal[[code]][,3],
             y=signal[[code]][,column],
             xlab='t(s)',ylab=label,pch=16,main=code)
    }
    tit <- paste0(names(x)[run],':',x[[run]]$name)
    mtext(tit, side=3, line=1.75, outer=TRUE)
    mtext(dat[[run]]$mdf$name, side=3, line=0.25, outer=TRUE)
}

plot.drift <- function(dc,run=1){
    meas <- dc[[run]]
    oldpar <- par(mfrow=c(2,2),mar=c(4,3.5,1,1),oma=c(0,0,3,0),mgp=c(2,1,0))
    on.exit(par(oldpar))
    FAR <- (meas$detectors %in% 'F')
    for (i in 1:ncol(meas$tim)){
        z <- ifelse(FAR[i],meas$ZF,meas$ZC)
        plot(x=meas$tim[,i],y=meas$sig[,i]-z,pch=16,
             xlab='t(s)',ylab=meas$labels[i])
        g <- meas$groups[i]
        b <- meas$alpha[i] + meas$gamma[g]*hours(meas$tim[,i])
        pred <- exp(b)/ifelse(FAR[i],1,meas$dwell[i])
        lines(x=meas$tim[,i],y=pred)
    }
    mtext(meas$name,side=3,line=1,outer=TRUE)
}
