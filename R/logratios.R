logratios <- function(bdat,num,den){
    err <- tab <- data.frame()
    for (i in 1:nrow(bdat$tab)){
        lab <- paste0('ln[',num,'/',den,']')
        tab[i,lab] <- bdat$tab[i,num] - bdat$tab[i,den]
        err[i,lab] <- bdat$err[i,num] - bdat$err[i,den]
    }
    tab$type <- bdat$type
    tab$comment <- bdat$comment
    list(tab=tab,err=err)
}
