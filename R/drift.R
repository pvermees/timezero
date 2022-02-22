zero <- function(adat){
    nms <- names(adat)
    nd <- length(adat)
    out <- list()
    for (nm in nms){
        meas <- adat[[nm]]
        FAR <- which(meas$detectors %in% 'F')
        SEM <- which(meas$detectors %in% 'C')
        sig <- meas$sig
        for (i in FAR){
            sig[,i] <- sig[,i] - meas$ZF
        }
        for (i in SEM){
            sig[,i] <- sig[,i] - meas$ZC
        }
        out[[nm]] <- list(tim=meas$tim,sig=sig,FAR=FAR,SEM=SEM)
    }
    out
}

drift <- function(dat){
    adat <- average.all.blocks(dat)
    out <- adat
    zdat <- zero(adat)
    nd <- length(adat)
    for (i in 1:nd){
        meas <- adat[[i]]
        ng <- length(unique(meas$groups))
        a <- log(zdat[[i]]$sig[1,])
        k <- log(colMeans(zdat[[i]]$sig))
        g <- rep(-4,ng)
        init <- c(a,k,g)
        fit <- optim(init,LL.drift,meas=meas)
        np <- length(fit$par)
        nak <- (np-ng)/2
        out[[i]]$alpha <- fit$par[1:nak]
        out[[i]]$kappa <- fit$par[(nak+1):(2*nak)]
        out[[i]]$gamma <- tail(fit$par,ng)
    }
    class(out) <- append('drift',class(out))
    out
}

LL.drift <- function(p,meas){
    FAR <- (meas$detectors %in% 'F')
    SEM <- (meas$detectors %in% 'C')
    ng <- length(unique(meas$groups))
    np <- length(p)
    nak <- (np-ng)/2
    alpha <- p[1:nak]
    kappa <- p[(nak+1):(2*nak)]
    gamma <- tail(p,n=ng)
    out <- 0
    for (i in 1:nak){
        g <- meas$groups[i]
        if (FAR[i]){
            pred <- meas$ZF + hall(alpha[i],gamma[g],kappa[i],meas,i)
            misfit <- pred - meas$sig[,i]
            out <- out + sum(misfit^2)/var(misfit)
        }
        if (SEM[i]){
            pred <- meas$ZC + hall(alpha[i],gamma[g],kappa[i],meas,i)
            obs <- round(meas$sig[,i] * meas$dwell[i])
            out <- out - sum(dpois(x=obs,lambda=pred,log=TRUE))
        }
    }
    out
}

drift.old <- function(dat){
    adat <- average.all.blocks(dat)
    out <- adat
    zdat <- zero(adat)
    nd <- length(adat)
    for (i in 1:nd){
        meas <- adat[[i]]
        ng <- length(unique(meas$groups))
        init <- c(log(zdat[[i]]$sig[1,]),rep(0,ng))
        fit <- optim(init,LL.drift,meas=meas)
        np <- length(fit$par)
        out[[i]]$alpha <- fit$par[1:(np-ng)]
        out[[i]]$gamma <- tail(fit$par,ng)
    }
    class(out) <- append('drift',class(out))
    out
}

LL.drift.old <- function(p,meas){
    FAR <- (meas$detectors %in% 'F')
    SEM <- (meas$detectors %in% 'C')
    ng <- length(unique(meas$groups))
    np <- length(p)
    alpha <- p[1:(np-ng)]
    gamma <- tail(p,n=ng)
    out <- 0
    for (i in 1:(np-ng)){
        g <- meas$groups[i]
        if (FAR[i]){
            pred <- meas$ZF + exp(alpha[i] + gamma[g]*hours(meas$tim[,i]))
            misfit <- pred - meas$sig[,i]
            out <- out + sum(misfit^2)/var(misfit)
        }
        if (SEM[i]){
            pred <- meas$ZC + exp(alpha[i] + gamma[g]*hours(meas$tim[,i]))
            obs <- round(meas$sig[,i] * meas$dwell[i])
            out <- out - sum(dpois(x=obs,lambda=pred,log=TRUE))
        }
    }
    out
}

hours <- function(tt){
    tt/3600
}
