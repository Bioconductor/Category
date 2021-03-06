##take an adjacency matrix, where the genes are columns, and a test
##stat vector and find the interesting sums
findAMstats = function(Amat, tstats) {
    pwLens = rowSums(Amat)
    expDE = Amat %*% tstats
    names(expDE) <- rownames(Amat)
    pw1s = pwLens > 1
    return(list(eDE = expDE[pw1s], lens=pwLens[pw1s]))
}

## Amat: adjacency matrix of gene -> categories
applyByCategory = function(stats, Amat, FUN=mean, ...)
{
  if(!(is.matrix(Amat)))
    stop("'Amat' needs to be a matrix.")
  if(ncol(Amat) != length(stats) )
    stop("'ncol(Amat)' needs to be the same as 'length(stats)'.")
  if(!is.logical(Amat)){
    if(!is.numeric(Amat))  
      stop("'Amat' needs to be logical or numeric.")
    Amat = (Amat!=0)
  }
  
  res = apply(Amat, 1, function(x) FUN(stats[x], ...))
  names(res) = rownames(Amat)
  return(res)
}

##given a set of AffyIDs, which have pathway data
probes2Path = function(pids, data="hgu133plus2") {
    pEnv = getAnnMap("PATH", data)
    inPW = mget(pids, pEnv, ifnotfound=NA)
    inPW[!is.na(inPW)]
}

getPathNames = function(iPW, organism = "hsa") {
    klist <- KEGGREST::keggList("pathway", organism = organism)
    knames <- paste0("path:", organism, iPW)
    kind <- lapply(stats::setNames(nm = knames), match, names(klist))
    kres <- vapply(kind, as.character, character(1L))
    kres[!is.na(kind)] <- klist[unlist(kind[!is.na(kind)])]
    lapply(strsplit(kres, " - "), `[`, 1L)
}

ttperm = function(x, fac, B=100, tsO=TRUE) {
    obs = rowttests(x, fac, tstatOnly= tsO)
    ans = vector(mode="list", length=B)
    for(i in 1:B) {
        p1 = sample(fac)
        ans[[i]] = rowttests(x, p1, tstatOnly = tsO)
    }
    return(list(obs=obs, perms=ans))
}

makeEBcontr = function(f1, hival) {
    if (!requireNamespace("EBarrays", quietly=TRUE))
      stop("need the EBarrays package for this feature")
    contr = rep(1, length(f1))
    p1 = paste(contr, collapse=",")
    contr2 = contr
    contr2[f1 == hival] = 2
    p2 = paste(contr2, collapse=",")
    EBarrays::ebPatterns(c(p1,p2))
}
