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

misfit.hall <- function(par,zdat,i,j){
    meas <- zdat[[i]]
    a <- par[1]
    g <- par[2]
    k <- par[3]
    tt <- hours(meas$tim[,j])
    obs <- meas$sig[,j]
    pred <- hall(a,g,k,tt)
    if (j %in% meas$FAR){
        out <- sum((obs - pred)^2)
    } else {
        O <- round((obs*meas$dwell[j] + meas$ZC*meas$zdwell))
        E <- pred*meas$dwell[j] + meas$ZC*meas$zdwell
        size <- sum(O)
        prob <- E/sum(E)
        out <- -sum(dpois(O,lambda=E,log=TRUE))
    }
    out
}
