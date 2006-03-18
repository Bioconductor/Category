##take an adjacency matrix, where the genes are columns, and a test
##stat vector and find the interesting sums
findAMstats = function(Amat, tstats) {
    pwLens = rowSums(Amat)
    expDE = Amat %*% tstats
    names(expDE) <- rownames(Amat)
    pw1s = pwLens > 1
    return(list(eDE = expDE[pw1s], lens=pwLens[pw1s]))
}

## Amat: either a 0/1 adjacency matrix of categories x genes,
##   or a graphNEL (directed) whose edges go from categories
##   to genes
applyByCategory = function(stats, Amat, FUN=mean, ...)
{
  if(is.matrix(Amat)) {
    if(ncol(Amat) != length(stats) )
        stop("wrong dimension for Amat")
    res = apply(Amat, 1, function(x)  FUN(stats[x!=0],...))
  } else {
    if(!inherits(Amat, "graphNEL"))
      stop("'Amat' must be a matrix or a graphNEL object.")
    if(is.null(names(stats)))
      stop("'stats' must be named, and names must correspond to the node names in 'Amat'")
    if(!identical(nodes(Amat), names(edges(Amat))))
      stop("'nodes(Amat)' and 'names(edges(Amat))' are not identical.")
    isGene = nodes(Amat) %in% names(stats)
    if(!any(isGene))
      warning("node names of 'Amat' do not match any of the names of 'stats'")
    res = sapply(edges(Amat)[!isGene], function(x) FUN(stats[x], ...))
  }
  return(res)  
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
