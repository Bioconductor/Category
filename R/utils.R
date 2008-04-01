getDataEnv <- function(name, lib) {
    lib <- sub(".db$", "", lib)
    get(paste(lib, name, sep=""))
}

getOrganism <- function(chip){
    paste("ORGANISM:", getAnnMap("ORGANISM", chip), sep="")
}
