readthedate <- function(x){
    thedate <- strptime(x,"%d-%b-%Y %H:%M")
    if (any(is.na(thedate))) thedate <- strptime(x,"%d-%m-%Y %H:%M")
    if (any(is.na(thedate))) thedate <- strptime(x,"%Y-%m-%d %H:%M:%S")
    if (any(is.na(thedate))) thedate <- strptime(x,"%d/%m/%y %H:%M")
    if (any(is.na(thedate))) thedate <- strptime(x,"%d-%b-%Y")
    return(as.numeric(thedate))
}
