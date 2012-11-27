getDataEnv <- function(name, lib) {
    get(paste(lib, name, sep=""))
}

## here getAnnMap is NOT returning a bimap.
## It's also not really accessing one.
getOrganism <- function(chip){
    paste("ORGANISM:", getAnnMap("ORGANISM", chip), sep="")
}
