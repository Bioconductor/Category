findLargest = function(gN, testStat, data="hgu133plus2") {
    LLe = get(paste(data, "LOCUSID", sep=""))
    lls = unlist(mget(gN, LLe))
    if(length(testStat) != length(gN) )
        stop("testStat and gN must be the same length")
    if( is.null(names(testStat)) )
        names(testStat) = gN
    tSsp = split.default(testStat, lls)
    sapply(tSsp, function(x) names(which.max(x)))
}

##take an adjacency matrix, where the genes are columns, and a test
##stat vector and find the interesting sums
findAMstats = function(Amat, tstats) {
    pwLens = rowSums(Amat)
    expDE = Amat %*% tstats
    pw1s = pwLens > 1
    return(list(eDE = expDE[pw1s], lens=pwLens[pw1s]))
}

##it seems rather simple to apply any function
applyByCategory = function(stats, Amat, FUN=mean, ...)
{
    if(ncol(Amat) != length(stats) )
        stop("wrong dimension for Amat")
    apply(Amat, 1, function(x)  FUN(stats[x==1],...))
}

##given a set of AffyIDs, which have pathway data
probes2Path = function(pids, data="hgu133plus2") {
    pEnv = get(paste(data, "PATH", sep=""))
    inPW = mget(pids, pEnv, ifnotfound=NA)
    inPW[!is.na(inPW)]
}

getPathNames = function(iPW)
    mget(iPW, KEGGPATHID2NAME, ifnotfound = NA)

ttperm = function(x, fac, B=100, tsO=TRUE) {
    obs = rowttests(x, fac, tstatOnly= tsO)
    ans=NULL
    for(i in 1:B) {
        p1 = sample(fac)
        ans[[i]] = rowttests(x, p1, tstatOnly = tsO)
    }
    return(list(obs=obs, perms=ans))
}

makeEBcontr = function(f1, hival) {
    contr = rep(1, length(f1))
    p1 = paste(contr, collapse=",")
    contr2 = contr; contr2[f1 == hival] = 2;
    p2 = paste(contr2, collapse=",")
    ebPatterns(c(p1,p2))
}
