hall <- function(a,g,k,tt){
    exp(a-g*tt) + (1-exp(-g*tt))*k/g
}

init.hall <- function(zdat,i){
    gmisfit <- function(g,a,tmiddle,cmiddle,tlast,clast){
        pred <- exp(a-g*tlast) +
            (cmiddle-exp(a-g*tmiddle))*(1-exp(-g*tlast))/(1-exp(-g*tmiddle))
        (clast - pred)^2
    }
    nr <- nrow(zdat$sig)
    pos <- zdat$sig[,i]>0
    first <- zdat$sig[pos,i][1]
    a <- log(first)
    clast <- zdat$sig[nr,i]
    tlast <- hours(zdat$tim[nr,i])
    cmiddle <- zdat$sig[round(nr/2),i]
    tmiddle <- hours(zdat$tim[round(nr/2),i])
    g <- optimise(gmisfit,interval=c(-10,10),
                  a,tmiddle,cmiddle,tlast,clast)$minimum
    k <- (cmiddle-exp(a-g*tmiddle))*g/(1-exp(-g*tmiddle))
    c(a,g,k)
}

misfit.hall <- function(par,zdat,j){
    n <- length(j)
    a <- par[1:n]
    g <- par[n+1]
    k <- par[(n+2):(2*n+1)]
    out <- 0
    for (i in 1:n){
        jj <- j[i]
        tt <- hours(zdat$tim[,jj])
        obs <- zdat$sig[,jj]
        pred <- hall(a[i],g,k[i],tt)
        if (jj %in% zdat$FAR){
            out <- sum((obs - pred)^2)
        } else {
            O <- round((obs*zdat$dwell[jj] + zdat$ZC*zdat$zdwell))
            E <- pred*zdat$dwell[jj] + zdat$ZC*zdat$zdwell
            size <- sum(O)
            prob <- E/sum(E)
            out <- out - sum(dpois(O,lambda=E,log=TRUE))
        }
    }
    out
}

hallfit <- function(zdat){
    groups <- unique(zdat$groups)
    out <- zdat
    out$agk <- list()
    class(out) <- 'hall'
    for (g in groups){
        j <- which(zdat$groups %in% g)
        n <- length(j)
        agk <- NULL
        for (jj in j){
            agk <- rbind(agk,init.hall(zdat,jj))
        }
        init <- c(agk[,1],agk[1,2],agk[,3])
        lower <- c(init[1:n]-1,init[n+1]-1,init[(n+2):(2*n+1)]/2)
        upper <- c(init[1:n]+1,init[n+1]+1,init[(n+2):(2*n+1)]*2)
        out$agk[[g]] <- optim(init,misfit.hall,method='L-BFGS-B',
                              lower=lower,upper=upper,zdat=zdat,j=j)$par
    }
    out
}

plot.hall <- function(zdat){
    n <- length(zdat$codes)
    nr <- ceiling(sqrt(n))
    nc <- ceiling(n/nr)
    oldpar <- par(mfrow=c(nr,nc))
    on.exit(par(oldpar))
    groups <- unique(zdat$groups)
    for (g in groups){
        agk <- zdat$agk[[g]]
        j <- which(zdat$groups %in% g)
        nj <- length(j)
        for (i in 1:nj){
            jj <- j[i]
            tt <- hours(zdat$tim[,jj])
            pred <- hall(agk[i],agk[nj+1],agk[nj+1+i],tt)
            obs <- zdat$sig[,jj]
            plot(tt,obs,ylab='y')
            lines(tt,pred)
        }
    }
}
