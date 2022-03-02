mdf <- function(fn){
    if (fn %in% c('He4He3_Small_Signals.mdf',
                  'He4He3_Small_Signals_Long.mdf',
                  'He4He3_Small_Signals_Very_Long.mdf')){
        out <- list(codes = c('101','102','103','104'),
                    labels = c('4He (V)','4He (cps)','3He (cps)','12C2+ (cps)'),
                    columns = c(1,2,2,2),
                    detectors = c('F','C','C','C'),
                    master = 1,
                    groups = c(1,1,1,2)
                    )
    } else if (fn %in% c('He4He3_Large_Signals_short.mdf',
                         'He4He3_Large_Signals.mdf',
                         'He4He3_Large_Signals_Long.mdf',
                         'He4He3_Large_Signals_Very_Long.mdf')){
        out <- list(codes = c('101','102','103'),
                    labels = c('4He (V)','3He (cps)','12C2+ (cps)'),
                    columns = c(1,2,2),
                    detectors = c('F','C','C'),
                    master = 1,
                    groups = c(1,1,2)
                    )
    }
    out
}

average.all.blocks <- function(dat){
    out <- list()
    for (nm in names(dat)){
        out[[nm]] <- average.blocks(dat[[nm]])
    }
    out
}

average.blocks <- function(meas){
    ZF <- mean(meas$signal[['1']][,1])
    ZC <- mean(meas$signal[['1']][,2])
    numblocks <- meas$mdf$numblocks
    method <- mdf(meas$mdf$name)
    codes <- method$codes
    numcodes <- length(codes)
    sig <- tim <- matrix(NA,numblocks,numcodes)
    colnames(sig) <- colnames(tim) <- codes
    for (code in codes){
        nperb <- nrow(meas$signal[[code]])/numblocks
        col <- method$columns[which(method$codes %in% code)]
        for (b in 1:numblocks){
            i <- (b-1)*nperb + (1:nperb)
            tim[b,code] <- mean(meas$signal[[code]][i,3])
            sig[b,code] <- mean(meas$signal[[code]][i,col])
        }
    }
    out <- method
    out$ZF <- ZF
    out$ZC <- ZC
    out$tim <- tim
    out$sig <- sig
    out$zdwell <- as.numeric(meas$mdf$zero[,5])
    out$dwell <- as.numeric(meas$mdf$cycle[,5])
    out$name <- meas$name
    out
}

naive.processing <- function(dat){
    adat <- average.all.blocks(dat)
    nd <- length(adat)
    out <- matrix(NA,nd,4)
    colnames(out) <- c('He3C/He4F','err','He3C/He4C','err')
    for (i in 1:nd){
        meas <- adat[[i]]
        if (ncol(meas$sig)<4){
            out[i,1:2] <- avg.ratio((meas$sig[,2]-meas$ZC)/
                                    (meas$sig[,1]-meas$ZF))
        } else {
            out[i,1:2] <- avg.ratio((meas$sig[,3]-meas$ZC)/
                                    (meas$sig[,1]-meas$ZF))
            out[i,3:4] <- avg.ratio((meas$sig[,3]-meas$ZC)/
                                    (meas$sig[,2]-meas$ZC))
        }
    }
    out
}

avg.ratio <- function(r){
    c(mean(r), sd(r)/sqrt(length(r)))
}
