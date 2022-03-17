get.hdat <- function(RUNdir,MDFdir,outliers){
    nms <- read.csv(paste0(RUNdir,'/names.csv'),as.is=TRUE)
    dat <- read.data(RUNdir,MDFdir,TRUE)
    adat <- average.all.blocks(dat)
    zdat <- zero(adat)
    hallall(zdat,outliers=outliers)
}
