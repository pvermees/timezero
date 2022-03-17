blank <- function(hdat,pairs){
    nms <- rownames(hdat$tab)
    isamps <- which(pairs$type %in% c(' s',' x'))
    err <- tab <- NULL
    lab <- colnames(hdat$tab)
    for (i in isamps){
        b1 <- pairs$blank1[i]
        b2 <- pairs$blank2[i]
        i1 <- which(nms %in% b1)
        if (b2 == ""){
            eb <- exp(hdat$tab[i1,])
        } else {
            i2 <- which(nms %in% b2)
            eb <- exp(hdat$tab[i1,])/2 + exp(hdat$tab[i2,])/2
        }
        tab <- rbind(tab, log( exp(hdat$tab[i,]) - eb ))
        err <- rbind(err, log( exp(hdat$tab[i,]) - eb ))
    }
    list(tab=as.data.frame(tab),
         err=as.data.frame(err),
         type=pairs$type[isamps],
         comment=pairs$comment[isamps])
}
