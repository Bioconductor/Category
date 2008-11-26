getDataEnv <- function(name, lib) {
    get(paste(lib, name, sep=""))
}

getOrganism <- function(chip){
    paste("ORGANISM:", getAnnMap("ORGANISM", chip), sep="")
}
