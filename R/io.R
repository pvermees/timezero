read.data <- function(RUNdir,MDFdir){ 
    out <- list()
    fnames <- Sys.glob(file.path(RUNdir,"*.RUN"))
    for (fn in fnames){
        con = file(fn, "r")
        meas <- list()
        meas$date <- read.date(con)
        meas$mdf <- read.mdf(con,MDFdir)
        line <- readLines(con,n=1) # the next line contains the sample name
        meas$name <- unlist(strsplit(line,'"'))[2]
        meas$signal <- read.signal(con)
        mname <- rev(strsplit(basename(fn),split="\\.")[[1]])[2]
        out[[mname]] <- meas
        close(con)
    }
    class(out) <- 'Noblesse'
    out
}

read.date <- function(con){
    while (TRUE){
        line <- readLines(con,n=1)
        if (length(line)==0) break
        if (grepl("Start Of Run",line)){
            txt <- unlist(strsplit(line,'#'))[2]
            return(readthedate(txt))
        }
    }
}

read.mdf <- function(con,MDFdir){
    while (TRUE){
        line <- readLines(con,n=1)
        if (length(line)==0) break
        if (grepl(".mdf",line,ignore.case=TRUE)) break
    }
    fname <- rev(strsplit(strsplit(line,"\"")[[1]][2],
                          "[\\]")[[1]])[1]
    out <- list(name=fname,zero=NULL,cycle=NULL)
    fn <- file.path(MDFdir,fname)
    MDFcon = file(fn, "r")
    block = readLines(MDFcon,n=13)
    nz <- as.numeric(strsplit(block[2],',')[[1]][2]) # number of zeros
    nc <- as.numeric(strsplit(block[3],',')[[1]][2]) # number of cycles
    block <- readLines(MDFcon,n=nz+nc)
    for (i in 1:(nz+nc)){
        txt <- strsplit(block[i],',')[[1]][-1]
        if (i>nz){
            txt[1] <- 100 + as.numeric(txt[1])
            out$cycle <- rbind(out$cycle,txt)
        } else {
            out$zero <- rbind(out$zero,txt)
        }
    }
    close(MDFcon)
    out
}

read.signal <- function(con){
    while (TRUE){
        line <- readLines(con,n=1)
        if (length(line)==0) break
        if (identical(line,'""')) break
    }
    out <- list()
    while (TRUE){
        line = readLines(con,n=1)
        if (length(line)==0) break
        txt <- unlist(strsplit(line,','))
        nt <- length(txt)
        hop <- txt[nt]
        if (!(hop %in% names(out))){
            out[[hop]] <- NULL
        }
        out[[hop]] <- rbind(out[[hop]],as.numeric(txt[-nt]))
    }
    out
}
