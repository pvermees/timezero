hall <- function(a,g,k,tt){
    exp(a-g*tt) + (1-exp(-g*tt))*k/g
}

init.hall <- function(zdat,i){
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

gmisfit <- function(g,a,tmiddle,cmiddle,tlast,clast){
    pred <- exp(a-g*tlast) +
        (cmiddle-exp(a-g*tmiddle))*(1-exp(-g*tlast))/(1-exp(-g*tmiddle))
    (clast - pred)^2
}
