mdf <- function(fn){
    if (fn %in% c('He4He3_Small_Signals.mdf',
                  'He4He3_Small_Signals_Long.mdf',
                  'He4He3_Small_Signals_Very_Long.mdf')){
        out <- list(codes = c('101','102','103','104'),
                    labels = c('4He (V)','4He (cps)','3He (cps)','12C2+ (cps)'),
                    columns = c(1,2,2,2),
                    detectors = c('F','C','C','C'),
                    master = 1,
                    groups = c(1,2,2,3)
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
                    groups = c(1,2,3)
                    )
    }
    out
}

naive.processing <- function(dat){
    x <- average.all.blocks(dat)
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
            sig[b,code] <- sum(meas$signal[[code]][i,col])
        }
    }
    out <- method
    out$ZF <- ZF
    out$ZC <- ZC
    out$tim <- tim
    out$sig <- sig
    out$dwell <- as.numeric(meas$mdf$cycle[,5])
    out
}
