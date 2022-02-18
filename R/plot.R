plot.Noblesse <- function(x,run=1){
    method <- mdf(dat[[run]]$mdf$name)
    signal <- x[[run]]$signal
    nm <- length(method$codes)
    nc <- ceiling(sqrt(nm))
    nr <- ceiling(nm/nc)
    oldpar <- par(mfrow=c(nr,nc),mar=c(3,4,1,1),oma=c(0,0,3,0))
    on.exit(par(oldpar))
    for (i in 1:nm){
        code <- method$codes[i]
        column <- method$columns[i]
        label <- method$labels[i]
        plot(x=signal[[code]][,3],
             y=signal[[code]][,column],
             xlab='t(s)',ylab=label,pch=16)
    }
    tit <- paste(names(x)[run],
                 '(',x[[run]]$name,
                 ') using ',
                 dat[[run]]$mdf$name)
    mtext(tit, side=3, line=1, outer=TRUE)
}
