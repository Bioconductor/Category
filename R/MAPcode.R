 probes2MAP <- function (pids, data = "hgu133plus2") {
    pEnv = get(paste(data, "MAP", sep = ""))
    inMAP = mget(pids, pEnv, ifnotfound = NA)
    inMAP[!is.na(inMAP)]
}

   MAPAmat <- function (data, minCount=5) {
    if (!is.character(data) || length(data) != 1)
        stop("wrong argument")
    dataE = as.list(get(paste(data, "MAP", sep = "")))
    LLs = getLL(names(dataE), data)
    goodLLs = !duplicated(LLs) & !is.na(LLs)
    use = dataE[goodLLs]
    names(use) = LLs[goodLLs]
    ##some are in more than one cytochrome band and we keep only the first
    use2 = sapply(use, function(x) x[1])
    cytoLL = split(names(use2), use2)
    cts = sapply(cytoLL, length)
    goodCts = cts >= minCount
    cytoLL = cytoLL[goodCts]
    uniqLL = unique(unlist(cytoLL))
    Amat = sapply(cytoLL, function(x) {
        mtch = match(x, uniqLL)
        zeros = rep(0, length(uniqLL))
        zeros[mtch] = 1
        zeros
    })
    dimnames(Amat) = list(uniqLL, names(cytoLL))
    return(Amat)
}

