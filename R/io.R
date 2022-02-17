read.data <- function(dname){ 
    out <- list()
    filenames <- Sys.glob(paste0(dname,"/*.RUN"))
    nf <- length(filenames)
    thedate <- rep(NA,nf)
    labels <- rep(NA,nf)
    irr <- rep(NA,nf)
    pos <- rep(NA,nf)
    hops=list('101'=c('He4[F]','b1'),
              '102'=c('b2','He4[C]'),
              '103'=c('b3','He3[C]'),
              '104'=c('b4','C12+++'))
            # '1'=c('Far','IC0')
    keep <- list('101'=c(TRUE,FALSE),
                 '102'=c(FALSE,TRUE),
                 '103'=c(FALSE,TRUE),
                 '104'=c(FALSE,TRUE))
    for (i in 1:nf){
        con = file(filenames[i], "r")
        while (TRUE){
            line = readLines(con,n=1)
            if (grepl("Start Of Run",line)){
                txt <- unlist(strsplit(line,'#'))[2]
                thedate[i] <- readthedate(txt)
            } else if (grepl(".mdf",line,ignore.case=TRUE)){
                line = readLines(con,n=1) # the next line contains the sample name
                sname <- unlist(strsplit(line,'"'))[2]
                naliquots <- length(which(grepl(sname,labels)))
                if (naliquots > 0) sname <- paste(sname,naliquots,sep=' ')
                labels[i] <- sname
            } else if (identical(line,'""')){ # start of signal
                out <- addSignal(out,con,hops,keep,thedate=thedate,
                                 thelabel=thelabel,theirr=NA,thepos=NA)
                break
            }
        }
        close(con)
    }
    for (hop in names(hops)){
        out[[hop]]$thedate <- thedate
        out[[hop]]$labels <- labels
        out[[hop]]$irr <- irr
        out[[hop]]$pos <- pos
        class(out[[hop]]) <- 'timeresolved'
    }
    class(out) <- 'Noblesse'
    out
}

addSignal <- function(dat,con,hops,keep,thedate=NA,
                      thelabel=NA,theirr=NA,thepos=NA){
    out <- dat
    thetime <- list()
    d <- list()
    while (TRUE){
        line = readLines(con,n=1)
        if (length(line)==0) break
        txt <- unlist(strsplit(line,','))
        nt <- length(txt)
        hop <- txt[nt]
        if (hop %in% names(hops)){
            thetime[[hop]] <- append(thetime[[hop]],as.numeric(txt[nt-1]))
            sweep <- txt[-nt]
            d[[hop]] <- rbind(d[[hop]],
                              as.numeric(sweep[keep[[hop]]]))
        }
    }
    for (hop in names(hops)){
        if (!(hop %in% names(out))){
            out[[hop]] <- list(masses=NULL,thetime=NULL,d=NULL)
        }        
        i <- keep[[hop]]
        out[[hop]]$masses <- hops[[hop]][i]
        out[[hop]]$thetime <- cbind(out[[hop]]$thetime,thetime[[hop]])
        out[[hop]]$d <- cbind(out[[hop]]$d,d[[hop]])
    }
    out
}
