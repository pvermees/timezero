hall <- function(a,g,k,tt){
    exp(a-g*tt) + (1-exp(-g*tt))*k/g
}

init.hall <- function(zdat,j,outliers=NULL){
    gmisfit <- function(g,a,tmiddle,cmiddle,tlast,clast){
        pred <- exp(a-g*tlast) +
            (cmiddle-exp(a-g*tmiddle))*(1-exp(-g*tlast))/(1-exp(-g*tmiddle))
        (clast - pred)^2
    }
    if (!is.null(outliers)){
        sig <- zdat$sig[-outliers,]
        tim <- zdat$tim[-outliers,]
    } else {
        sig <- zdat$sig
        tim <- zdat$tim
    }
    nr <- nrow(sig)
    pos <- sig[,j]>0
    first <- sig[pos,j][1]
    a <- log(first)
    clast <- sig[nr,j]
    tlast <- hours(tim[nr,j])
    cmiddle <- sig[round(nr/2),j]
    tmiddle <- hours(tim[round(nr/2),j])
    g <- optimise(gmisfit,interval=c(-10,10),
                  a,tmiddle,cmiddle,tlast,clast)$minimum
    k <- (cmiddle-exp(a-g*tmiddle))*g/(1-exp(-g*tmiddle))
    c(a,g,k)
}

misfit.hall <- function(par,zdat,j,outliers=NULL){
    n <- length(j)
    a <- par[1:n]
    g <- par[n+1]
    k <- par[(n+2):(2*n+1)]
    out <- 0
    if (!is.null(outliers)){
        sig <- zdat$sig[-outliers,]
        tim <- zdat$tim[-outliers,]
    } else {
        sig <- zdat$sig
        tim <- zdat$tim
    }
    for (i in 1:n){
        jj <- j[i]
        tt <- hours(tim[,jj])
        obs <- sig[,jj]
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

hallfit <- function(zdat,outliers=NULL){
    out <- zdat
    class(out) <- 'hall'
    out$agk <- list()
    n <- length(zdat$codes)
    for (j in 1:n){
        init <- init.hall(zdat,j,outliers=outliers)
        lower <- c(init[1]-1,init[2]-1,init[3]/2)
        upper <- c(init[1]+1,init[2]+1,init[3]*2)
        code <- zdat$codes[j]
        out$agk[[code]] <- optim(init,misfit.hall,method='L-BFGS-B',
                                 lower=lower,upper=upper,
                                 zdat=zdat,j=j,outliers=outliers)$par
    }
    out
}
plot.hall <- function(zdat){
    n <- length(zdat$codes)
    nr <- ceiling(sqrt(n))
    nc <- ceiling(n/nr)
    oldpar <- par(mfrow=c(nr,nc))
    on.exit(par(oldpar))
    for (j in 1:n){
        agk <- zdat$agk[[j]]
        tt <- hours(zdat$tim[,j])
        pred <- hall(agk[1],agk[2],agk[3],tt)
        obs <- zdat$sig[,j]
        plot(tt,obs,ylab='y')
        lines(tt,pred)
    }
}

hallfitjoint <- function(zdat){
    groups <- unique(zdat$groups)
    out <- zdat
    out$agk <- list()
    class(out) <- 'hall2'
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
plot.halljoint <- function(zdat){
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
