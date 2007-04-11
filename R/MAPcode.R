## Since we now depends on the Matrix package, and Matrix
## does evil things to base::cbind, we make our own copy
## for basic uses.  The cbind that results from loading Matrix
## is unusable for perf reasons.
base_cbind <- function (..., deparse.level = 1) 
  .Internal(cbind(deparse.level, ...))

probes2MAP <- function (pids, data = "hgu133plus2") {
    pEnv = get(paste(data, "MAP", sep = ""))
    inMAP = mget(pids, pEnv, ifnotfound = NA)
    inMAP[!is.na(inMAP)]
}

   MAPAmat <- function (data, minCount=5) {
    if (!is.character(data) || length(data) != 1)
        stop("wrong argument")
    dataE = as.list(getAnnMap("MAP", chip=data))
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


### ------------------------------------------------------------------

cleanMapAndsOrs <- function(p2m) {
    ## The map data is a bit odd and we get things like: "Xq28 and Yq12"
    ## and "Xp22.3 or Yp11.3".  This function maps these to a proper
    ## vector of length 2.
    lapply(p2m, function(cbs) {
        unlist(strsplit(cbs, " and | or "))
    })
}


cleanMapWeird <- function(p2m) {
    ## Current MAP annotation data contains weird entries like:
    ##     "1 69.9 cM"
    ##     "1 E4"
    ##     "11 29.0 cM"
    ##     "11 B1.3"
    ##     "19 23.0 cM"
    ##     "19 C1"
    ##     "1q23.1 according to Sierra (Genomics 79, 177, 2002) [AFS]"
    ##     "2 A1"
    ##
    ## This function removes these for now.  The strategy is simply to
    ## look for annotations with spaces and drop them.
    lapply(p2m, function(cbs) {
        badIdx <- grep(" ", cbs)
        if (length(badIdx) > 0) {
            good <- cbs[-badIdx]
            if (length(good) == 0)
              good <- as.character(NA)
        } else {
            good <- cbs
        }
        good
    })
}


cleanRanges <- function(probe2chr) {
    ## We use (p|q) even though it will match nonsense annotations like
    ## "6p23-q24".  These seem to appear in the data.  But it is OK
    ## since we are taking the longest common suffix and that will just
    ## be the chromosome in these cases.
    rangePat <- "[0-9XY]+(q|p)([0-9.]+|ter)-(q|p)([0-9.]+|ter)"
    cenPat <- "cen"

    probe2chr <- lapply(probe2chr, function(x) {
        ## this could be made faster
        sapply(x, function(z) {
            if (length(grep("-", z, fixed=TRUE))) {
                if (length(grep(rangePat, z, perl=TRUE))) {
                    arm.loc <- gregexpr("(q|p|-)", z)[[1L]]
                    lhs <- substr(z, 1L, arm.loc[2L]-1L)
                    rhs <- paste(substr(z, 1L, arm.loc[1L]-1L),
                                 substr(z, arm.loc[2L]+1L, nchar(z)), sep="")
                    ans <- lcPrefixC(c(lhs, rhs))
                    nc <- nchar(ans)
                    if (substr(ans, nc, nc) == ".")
                      ans <- substr(ans, 1L, nc - 1L)
                    ans
                } else if (length(grep(cenPat, z, perl=TRUE))) {
                    ## FIXME:
                    ## for range including centromere, we just
                    ## take the entire chromosome.
                    ## Could we do better and include arm info?
                    substr(z, 1L, regexpr("cen|q|p", z) - 1L)
                } else {
                    parts <- strsplit(z, "-", fixed=TRUE)[[1]]
                    ans <- lcPrefix(parts)
                    if (nchar(ans) < 1)
                      ans <- substr(z, 1L, regexpr("cen|q|p", z) - 1L)
                    warning("unexpected band notation: ", z, " using ", ans,
                            call.=FALSE)
                    ans
                }
            } else {
                z
            }
        })
    })
    probe2chr
}


.chrMapToEntrez <- function(m2p, entrezMap, univ) {
    allEntrez <- unique(unlist(as.list(entrezMap)))
    allEntrez <- allEntrez[!is.na(allEntrez)]
    if (!is.null(univ) && length(univ) > 0)
      univ <- intersect(univ, allEntrez)
    else
      univ <- allEntrez
    m2eg <- lapply(m2p, function(x) {
        x <- unique(x[!is.na(x)])
        if (length(x) == 0)
          eg <- NULL
        else {
            eg <- unlist(mget(x, entrezMap, ifnotfound=NA))
            eg <- intersect(unique(eg[!is.na(eg)]), univ)
            if (length(eg) == 0)
              eg <- NULL
        }
        eg
    })
    notNull <- sapply(m2eg, function(x) !is.null(x))
    m2eg[notNull]
}


annBaseName <- function(p) {
    baseName <- p@datPkg@baseName
}

makeChrMapToEntrez <- function(chip, univ) {
    .getMap <- function(map)
      getAnnMap(map=map, chip=chip)
    aData <- .getMap(map="MAP")
    probe2chr <- as.list(aData)
    ## remove NAs
    probeNA <- sapply(probe2chr, function(x) {
        length(x) == 1 && is.na(x)
    })
    probe2chr <- probe2chr[!probeNA]
    probe2chr <- cleanMapAndsOrs(probe2chr)
    probe2chr <- cleanMapWeird(probe2chr)
    probe2chr <- cleanRanges(probe2chr)

    lens <- listLen(probe2chr)
    chr <- unlist(probe2chr)
    pbs <- rep(names(probe2chr), lens)
    m2p <- split(pbs, chr)
    ## XXX: remove all annotations that remain that have '-'
    hasRange <- grep("-", names(m2p))
    if (length(hasRange)) {
        warning("dropping ", length(hasRange), " items with weird range")
        m2p <- m2p[-hasRange]
    }
    onlyDot <- grep("^\\.$", names(m2p))
    if (length(onlyDot))
      m2p <- m2p[-onlyDot]
    invalid <- grep("[^0-9.qpXYter]+", names(m2p))
    if (length(invalid)) {
        warning("dropping invalid items: ",
                paste(names(m2p)[invalid], collapse=", "))
        m2p <- m2p[-invalid]
    }
    .chrMapToEntrez(m2p, .getMap("ENTREZID"), univ)
}


parseChrMap <- function(x) {
    ## Given a chromosome band annotation (see examples below),
    ## return a vector giving the path from the root:
    ##
    ## 20p12.2 => 20, 20p, 20p1, 20p12, 20p12.2
    ##
    ## x1 <- "2q32"
    ## x2 <- "4q21.22"
    ## x3 <- "20p12.2"
    ## x4 <- "2qter"
    ##
    arm.pos <- regexpr("p|q", x)
    if (arm.pos < 0)
      return(x)
    chr <- substr(x, 1, arm.pos - 1)
    chrArm <- substr(x, 1, arm.pos)
    if (substr(x, arm.pos + 1, arm.pos + 3) == "ter")
      return(c(chr, chrArm, x))
    if (nchar(x) == nchar(chrArm))
      return(c(chr, chrArm))
    sbs <- strsplit(substr(x, arm.pos+1, nchar(x)), "")[[1]]
    bands <- character(length(sbs)+1)
    bands[1] <- chrArm
    i <- 1
    j <- 2
    prev <- bands[1]
    while (TRUE) {
        if (sbs[i] == ".") {
            prev <- paste(prev, ".", sep="")
            i <- i + 1
            next
        }
        bands[j] <- paste(prev, sbs[i], sep="")
        prev <- bands[j]
        i <- i + 1
        j <- j + 1
        if (i > length(sbs))
          break
    }
    if (i == j)                         # there was a '.'
      c(chr, bands[1:(i-1)])
    else
      c(chr, bands)
}


makeChrBandGraph <- function(chip, univ=NULL) {
    m2p <- makeChrMapToEntrez(chip, univ)
    allMaps <- lapply(names(m2p), parseChrMap)
    vv <- lapply(allMaps, function(x) {
        L <- length(x)
        ## XXX: some will have L == 1, will induce NA's
        rbind(x[1:(L - 1)], x[2:L])
    })
    vvAll <- do.call(base_cbind, vv)
    vvStr <- paste(vvAll[1, ], vvAll[2, ], sep="+")
    ## remove duplicate edges
    dupIdx <- which(duplicated(vvStr))
    if (length(dupIdx) > 0 && any(dupIdx > 0))
      vvUniq <- vvAll[, -dupIdx]
    else
      vvUniq <- vvAll
    dimnames(vvUniq) <- list(c("from", "to"), NULL)
    vvT <- t(vvUniq)
    allNodes <- unique(unlist(allMaps))
    ## remove NA's
    NAitems <- which(is.na(vvT[, 2]))
    if (length(NAitems) > 0 && any(NAitems > 0))
      vvT <- vvT[-NAitems, ]
    ## remove self-loops
    selfLoops <- which(vvT[, 1] == vvT[, 2])
    if (length(selfLoops) > 0 && any(selfLoops > 0))
      vvT <- vvT[-selfLoops, ]
    ## add root node
    org <- paste("ORGANISM:", getAnnMap("ORGANISM", chip), sep="")
    chrNames <- names(getAnnMap("CHRLENGTHS", chip))
    orgLinks <- base_cbind(rep(org, length(chrNames)), chrNames)
    vvT <- rbind(orgLinks, vvT)
    g <- ftM2graphNEL(vvT, edgemode="directed")
    g <- addChrBandAnnotation(g, m2p)
    g
}

addChrBandAnnotation <- function(g, m2p) {
    nodeDataDefaults(g, "geneIds") <- as.character(NA)
    bands <- rev(tsort(g)[-1])          # remove root node
    gEdges <- edges(g)
    for (n in bands) {
        ids <- m2p[[n]]
        kids <- gEdges[[n]]
        if (length(kids)) {
            ids <- unique(c(ids, unlist(nodeData(g, n=kids, "geneIds"))))
        }
        if (length(ids) > 0)
          nodeData(g, n, "geneIds") <- list(ids[!is.na(ids)])
    }
    g
}

makeChrBandInciMat <- function(chrGraph) {
    gnodes <- nodes(chrGraph)
    rootIdx <- grep("^ORGANISM", gnodes)
    v1 <- nodeData(chrGraph, n=gnodes[-rootIdx] , attr="geneIds")
    v1 <- v1[sapply(v1, function(z) !any(is.na(z)))]
    v1 <- lapply(v1, as.character)
    allGeneIDs <- unique(unlist(v1))
    nr <- length(allGeneIDs)
    nc <- length(v1)
    mat <- matrix(0L, nrow=nr, ncol=nc)
    dimnames(mat) <- list(allGeneIDs, names(v1))
    for (cb in names(v1)) {
        mat[v1[[cb]], cb] <- 1L         # XXX: beware partial match!
    }
    mat
}

chrBandInciMat <- function(chip, univ=NULL) {
    makeChrBandInciMat(makeChrBandGraph(chip=chip, univ=univ))
}
