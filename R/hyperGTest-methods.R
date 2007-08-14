setMethod("hyperGTest",
          signature(p="HyperGParams"),
          function(p) .hyperGTestInternal(p))

.hyperGTestInternal <- function(p, className="HyperGResult") {
    p <- makeValidParams(p)
    origGeneIds <- geneIds(p)
    universeGeneIds(p) <- universeBuilder(p)
    selected <- intersect(geneIds(p), universeGeneIds(p))
    geneIds(p) <- selected
    cat2Entrez <- categoryToEntrezBuilder(p)
    stats <- .doHyperGTest(p, cat2Entrez, list(),
                           selected)
    ord <- order(stats$p)
    new(className,
        pvalues=stats$p[ord],
        oddsRatios=stats$odds[ord],
        expectedCounts=stats$expected[ord],
        catToGeneId=cat2Entrez[ord],
        annotation=annotation(p),
        geneIds=geneIds(p),
        testName=categoryName(p),
        pvalueCutoff=pvalueCutoff(p),
        testDirection=testDirection(p))
}

chrMap_hg_test <- function(p) {
    chrGraph <- p@chrGraph

    nodeDataDefaults(chrGraph, "pvalue") <- 1
    nodeDataDefaults(chrGraph, "geneIds") <- numeric(0)
    nodeDataDefaults(chrGraph, "condGeneIds") <- numeric(0)
    nodeDataDefaults(chrGraph, "oddsRatio") <- numeric(0)
    nodeDataDefaults(chrGraph, "expCount") <- numeric(0)

    ## now iterate leaves first doing tests and conditioning
    ## on all significant children.
    ## FIXME: consider replacing with RBGL tsort?
    needsProc <- chrGraph
    complete <- character(0)
    SIGNIF <- p@pvalueCutoff
    cat2Entrez <- nodeData(chrGraph, attr="geneIds")
    while (length(nodes(needsProc))) {
        numKids <- sapply(edges(needsProc), length)
        noKids <- names(numKids[numKids == 0])
        curCat2Entrez <- cat2Entrez[noKids]
        if (p@conditional) {
            curCatKids <- edges(chrGraph)[names(curCat2Entrez)]
            curCatKids <- removeLengthZero(curCatKids)
            if (length(curCatKids)) { ## sanity check
                ## they should be all complete
                stopifnot(all(unlist(curCatKids) %in% complete))
            }
            curCat2Entrez <- removeSigKidGenes(curCatKids, chrGraph,
                                               curCat2Entrez,
                                               SIGNIF, cat2Entrez)
            ## Store the conditioned cat => entrez map
            nodeData(chrGraph, n=names(curCat2Entrez),
                     attr="condGeneIds") <- curCat2Entrez
        }
        stats <- .doHyperGTest(p, curCat2Entrez, cat2Entrez,
                               p@geneIds)

        ## store the pvals, mark these nodes as complete,
        ## then compute the next set of nodes to do.
        noKids <- names(curCat2Entrez)
        ## drop names on pvals to avoid weird names upon unlisting
        nodeData(chrGraph, n=noKids,
                 attr="pvalue") <- as.numeric(stats$p)
        nodeData(chrGraph, n=noKids,
                 attr="oddsRatio") <- as.numeric(stats$odds)
        nodeData(chrGraph, n=noKids,
                 attr="expCount") <- as.numeric(stats$expected)
        complete <- c(complete, noKids)
        hasKids <- names(numKids[numKids > 0])
        needsProc <- subGraph(hasKids, needsProc)
    } ## end while
    p@chrGraph <- chrGraph
    p
}


removeLengthZero <- function(x) {
    wanted <- sapply(x, function(z) length(z) > 0)
    x[wanted]
}

removeSigKidGenes <- function(curCatKids, goDag, curCat2Entrez, SIGNIF,
                              cat2Entrez) {
    if (length(curCatKids)) {
        ## keep only those kids with SIGNIF pvalue
        curCatKids <- lapply(curCatKids, function(x) {
            pvKids <- nodeData(goDag, n=x, attr="pvalue")
            idx <- which(pvKids < SIGNIF)
            if (length(idx))
              x[idx]
            else
              character(0)
        })
        curCat2EntrezCond <- list()
        for (goid in names(curCat2Entrez)) {
            ## remove entrez ids that came from
            ## SIGNIF children
            kids <- curCatKids[[goid]]
            if (length(kids)) {
                kidEgs <- unlist(cat2Entrez[kids])
                newEgs <- setdiff(curCat2Entrez[[goid]], kidEgs)
                ## newEgs may be length 0
                curCat2EntrezCond[[goid]] <- newEgs
            } else {
                curCat2EntrezCond[[goid]] <- curCat2Entrez[[goid]]
            }
        }
        curCat2Entrez <- curCat2EntrezCond
    }
    curCat2Entrez
}



setMethod("hyperGTest",
          signature(p="ChrMapHyperGParams"),
          function(p) {
              p <- makeValidParams(p)
              if (numNodes(p@chrGraph) == 0)
                p@chrGraph <- makeChrBandGraph(annotation(p), p@universeGeneIds)

              univ <- unlist(nodeData(p@chrGraph, attr="geneIds"))
              univ <- unique(univ)
              p@universeGeneIds <- univ
              ## preserve names on geneIds
              p@geneIds <- p@geneIds[p@geneIds %in% p@universeGeneIds]

              p <- chrMap_hg_test(p)
              pvals <- unlist(nodeData(p@chrGraph, attr="pvalue"))
              pvord <- order(pvals)
              new("ChrMapHyperGResult",
                  chrGraph=p@chrGraph,
                  annotation=p@annotation,
                  geneIds=p@geneIds,
                  testName=categoryName(p),
                  testDirection=p@testDirection,
                  pvalueCutoff=p@pvalueCutoff,
                  conditional=p@conditional,
                  pvalue.order=pvord)
          })


geneGoHyperGeoTest <- function(entrezGeneIds, lib, ontology, universe=NULL)
{
    .Defunct("hyperGTest")
    if (missing(universe) || is.null(universe))
      universe <- character(0)
    params <- new("GOHyperGParams",
                  geneIds=entrezGeneIds,
                  universeGeneIds=universe,
                  annotation=lib,
                  ontology=ontology)
    hyperGTest(params)
}


geneKeggHyperGeoTest <- function(entrezGeneIds, lib, universe=NULL)
{
    .Defunct("hyperGTest")
    if (missing(universe) || is.null(universe))
      universe <- character(0)
    params <- new("KEGGHyperGParams",
                  geneIds=entrezGeneIds,
                  universeGeneIds=universe,
                  annotation=lib)
    hyperGTest(params)
}

.doHyperGInternal <- function(numW, numB, numDrawn, numWdrawn, over) {
    n21 <- numW - numWdrawn
    n12 <- numDrawn - numWdrawn
    n22 <- numB - n12

    odds_ratio <-  (numWdrawn * n22) / (n12 * n21)

    expected <- (numWdrawn + n12) * (numWdrawn + n21)
    expected <- expected / (numWdrawn + n12 + n21 + n22)

    if (over) {
        ## take the -1 because we want evidence for as extreme or more
        pvals <- phyper(numWdrawn - 1L, numW, numB,
                        numDrawn, lower.tail=FALSE)
    } else {
        pvals <- phyper(numWdrawn, numW, numB,
                        numDrawn, lower.tail=TRUE)
    }
    list(p=pvals, odds=odds_ratio, expected=expected)
}


.doHyperGTest <- function(p, curCat2Entrez, cat2Entrez, selected) {
    ## Here is how we conceptualize the test:
    ##
    ## The urn contains genes from the gene universe.  Genes annotated at a
    ## given cateogry term are white and the rest black.
    ##
    ## The number drawn is the size of the selected gene list.  The
    ## number of white drawn is the size of the intersection of the
    ## selected list and the genes annotated at the category.
    ##
    ## In the conditional case, currently only implemented for GO, the
    ## category ID annotation set has been reduced and we also adjust the
    ## selected list (num drawn) and gene universe according to what was
    ## removed by the conditioning.
    ##
    ## Here's a diagram based on using GO as the category:
    ##
    ##          inGO    notGO
    ##          White   Black
    ## selected  n11     n12
    ## not       n21     n22
    ##
    if (conditional(p)) {
        cat2RemovedEntrez <- lapply(names(curCat2Entrez),
                                    function(goid) {
                                        setdiff(cat2Entrez[[goid]],
                                                curCat2Entrez[[goid]])
                                    })

        ## White balls removed from urn by conditioning
        numSelectedRemoved <- sapply(cat2RemovedEntrez,
                                     function(x) sum(selected %in% x))

        numUnivRemoved <- listLen(cat2RemovedEntrez)
        ## Black balls removed from urn by conditioning
        numOtherRemoved <- numUnivRemoved - numSelectedRemoved
    } else {
        numSelectedRemoved <- rep(0, length(curCat2Entrez))
        numOtherRemoved <- numUnivRemoved <- numSelectedRemoved
    }

    ## Num white drawn (n11)
    numWdrawn <- sapply(curCat2Entrez,
                        function(x) sum(selected %in% x))

    ## Num white
    numW <- listLen(curCat2Entrez)

    ## Num black
    numB <- length(universeGeneIds(p)) - numUnivRemoved - numW

    ## Num drawn
    numDrawn <- length(selected) - numSelectedRemoved

    over <- testDirection(p) == "over"
    .doHyperGInternal(numW, numB, numDrawn, numWdrawn, over)
}
