getDataEnv <- function(name, lib) {
    get(paste(lib, name, sep=""), mode="environment")
}
