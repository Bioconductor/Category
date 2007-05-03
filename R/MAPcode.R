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

cleanMapAndsOrs <- function(e2m) {
    ## The map data is a bit odd and we get things like: "Xq28 and Yq12"
    ## and "Xp22.3 or Yp11.3".  This function maps these to a proper
    ## vector of length 2.
    for (eg in ls(e2m)) {
        e2m[[eg]] <- unlist(strsplit(e2m[[eg]], " and | or "))
    }
    e2m
}


cleanMapWeird <- function(e2m) {
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
    for (eg in ls(e2m)) {
        val <- e2m[[eg]]
        badIdx <- grep(" ", val)
        if (length(badIdx) > 0) {
            good <- val[-badIdx]
            if (length(good) == 0)
              good <- as.character(NA)
        } else {
            good <- val
        }
        e2m[[eg]] <- good
    }
    e2m
}


cleanRanges <- function(e2m) {
    ## We use (p|q) even though it will match nonsense annotations like
    ## "6p23-q24".  These seem to appear in the data.  But it is OK
    ## since we are taking the longest common suffix and that will just
    ## be the chromosome in these cases.
    rangePat <- "[0-9XY]+(q|p)([0-9.]+|ter)-(q|p)([0-9.]+|ter)"
    cenPat <- "cen"

    for (eg in ls(e2m)) {
        e2m[[eg]] <- sapply(e2m[[eg]], function(z) {
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
    }
    e2m
}


annBaseName <- function(p) {
    baseName <- p@datPkg@baseName
}

makeChrMapToEntrez <- function(chip, univ) {
    .getMap <- function(map)
      getAnnMap(map=map, chip=chip)
    probe2chr <- .getMap(map="MAP")
    if (!is.environment(probe2chr))
      probe2chr <- l2e(as.list(probe2chr))
    egs <- unique(unlist(mget(ls(probe2chr), .getMap(map="ENTREZID"))))
    egs <- egs[!is.na(egs)]
    if (!is.null(univ))
      egs <- intersect(egs, univ)
    egs <- as.character(egs)
    eg2chr <- new.env(parent=emptyenv(), hash=TRUE,
                      size=as.integer(1.20 * length(egs)))
    ## XXX: need to define a revmap method for environments
    eg2p <- l2e(mget(egs, revmap(.getMap(map="ENTREZID"))))
    for (eg in egs) {
        bands <- mget(eg2p[[eg]], probe2chr)
        bands <- bands[!is.na(bands)]
        if (length(bands))
          eg2chr[[eg]] <- unique(unlist(bands))
    }
    eg2chr <- cleanMapAndsOrs(eg2chr)
    eg2chr <- cleanMapWeird(eg2chr)
    eg2chr <- cleanRanges(eg2chr)

    m2eg <- reverseSplit(as.list(eg2chr))
    
    ## XXX: remove all annotations that remain that have '-'
    hasRange <- grep("-", names(m2eg))
    if (length(hasRange)) {
        warning("dropping ", length(hasRange), " items with weird range")
        m2eg <- m2eg[-hasRange]
    }
    onlyDot <- grep("^\\.$", names(m2eg))
    if (length(onlyDot))
      m2eg <- m2eg[-onlyDot]
    invalid <- grep("[^0-9.qpXYter]+", names(m2eg))
    if (length(invalid)) {
        warning("dropping invalid items: ",
                paste(names(m2eg)[invalid], collapse=", "))
        m2eg <- m2eg[-invalid]
    }
    m2eg
}


cb_parse_band_hsa <- function(x) {
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
    org <- paste("ORGANISM:", getAnnMap("ORGANISM", chip), sep="")
    if (org != "ORGANISM:Homo sapiens")
      stop("makeChrBandGraph can only deal with 'Homo sapiens' annotation",
           " found ", sQuote(org))
    m2p <- makeChrMapToEntrez(chip, univ)
    allMaps <- lapply(names(m2p), cb_parse_band_hsa)
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
    ## TODO, FIXME: should the return be transposed?
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

ann_list_to_mat <- function(L) {
    gsets <- names(L)
    nr <- length(gsets)
    allgenes <- unique(unlist(L))
    nc <- length(allgenes)
    M <- matrix(0L, nrow=nr, ncol=nc,
                dimnames=list(gsets, allgenes))
    for (gs in gsets) {
        M[gs, L[[gs]]] <- 1L
    }
    M
}

cb_childAnnList <- function(n, g) {
    ## n - node label, e.g., 12p1
    ## g - graph object representing the tree of chrom bands
    ##
    ## Returns a list named by children of n, values
    ## are vectors of gene IDs.  If n is not a node in
    ## g or has no children, return list()
    if (!(n %in% nodes(g)))
      return(list())
    kids <- edges(g)[[n]]
    if (length(kids))
      nodeData(g, n=kids, attr="geneIds")
    else
      list()
}

cb_childInciMat <- function(n, g) {
    ## Return an incidence matrix for the gene sets
    ## (rows) which are the children of node n in graph
    ## g.  The columns are the gene IDs
    dat <- cb_childAnnList(n, g)
    ann_list_to_mat(dat)
}

cb_childContinTable <- function(imat, selids, min.found=1L, min.k=1L) {
    ## Return a 2 x k contingency table
    ## First row is selected, second is not.  Columns
    ## are gene sets (rows of imat).
    ## If no gene set has greater or equal than min.found
    ## in the selids, then return integer(0).
    ##
    if (any(dim(imat) == 0))
      return(integer(0))
    if (nrow(imat) < min.k)
      return(integer(0))
    wh <- match(selids, colnames(imat), 0)
    submat <- imat[ , wh, drop=FALSE]
    if (ncol(submat) < min.found)       # no chance
      return(integer(0))
    tot <- rowSums(imat)
    found <- rowSums(submat)
    if (max(found) < min.found)
      return(integer(0))
    tab <- t(cbind(found, tot - found))
    dimnames(tab) <- list(c("selected", "not"), rownames(imat))
    tab
}

cb_contingency <- function(selids, chrList, chrGraph, testFun=chisq.test,
                           min.found=5L, min.k=1L) {
    chrMats <- lapply(chrList, cb_childInciMat, chrGraph)
    conTables <- lapply(chrMats, cb_childContinTable,
                        selids, min.found=min.found, min.k=min.k)
    conTables <- conTables[listLen(conTables) > 0]
    lapply(conTables, function(tab) {
        list(table=tab, result=testFun(tab))
    })
}

cb_sigBands <- function(b, p.value=0.01) {
    bands <- lapply(b, function(x) {
        if (x$result$p.value < p.value)
          colnames(x$table)
        else
          as.character(NA)
    })
    ans <- unlist(bands[!is.na(bands)])
    if (!length(ans))
      character(0)
    else
      ans
}

## ----

## This is another non-incidence matrix approach to
## generating the contingency tables.  Probably can
## junk this.
makeContinTable <- function(selids, ann_dat) {
    ## selids - vector of Gene IDs
    ## ann_dat - a list of vectors of Gene IDs; each
    ##           represents the IDs in a cateogry.
    ##           We assume the categories are close
    ##           to a partition of the genes.
    ##
    ## Returns a matrix with rows 'selected' and 'not' and a column for
    ## each element of ann_dat giving the counts for the contingency
    ## table.
    do.call(cbind, lapply(ann_dat, function(x) {
        found <- sum(selids %in% x)
        c(selected=found, not=(length(x) - found))
    }))
}
cb_makeContinTable <- function(n, g, selids) {
    dat <- cb_childAnnList(n, g)
    if (length(dat))
      makeContinTable(selids, dat)
    else
      dat
}

## ----
