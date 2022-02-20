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
