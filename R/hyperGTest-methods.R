setMethod("hyperGTest",
          signature(p="HyperGParams"), 
          function(p) {
              ##FIXME: add code for over/under representation handling
              ## also, reorg p-value calculation as in GOstats
              if (testDirection(p) != "over")
                stop("unsupported test direction: ", testDirection(p))
              origGeneIds <- geneIds(p)
              universeGeneIds(p) <- universeBuilder(p)
              selected <- intersect(geneIds(p), universeGeneIds(p))
              geneIds(p) <- selected
              cat2Entrez <- categoryToEntrezBuilder(p)
              stats <- .doHyperGTest(p, cat2Entrez, list(),
                                     selected)
              ord <- order(stats$p)
              new("HyperGResult",
                  pvalues=stats$p[ord],
                  oddsRatios=stats$odds[ord],
                  expectedCounts=stats$expected[ord],
                  catToGeneId=cat2Entrez[ord],
                  annotation=annotation(p),
                  geneIds=geneIds(p),
                  testName=categoryName(p),
                  pvalueCutoff=pvalueCutoff(p),
                  testDirection=testDirection(p))
          })


geneGoHyperGeoTest <- function(entrezGeneIds, lib, ontology, universe=NULL)
{
    .Deprecated("hyperGTest")
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
    .Deprecated("hyperGTest")
    if (missing(universe) || is.null(universe))
      universe <- character(0)
    params <- new("KEGGHyperGParams",
                  geneIds=entrezGeneIds,
                  universeGeneIds=universe,
                  annotation=lib)
    hyperGTest(params)
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
    if (isConditional(p)) {
        cat2RemovedEntrez <- lapply(names(curCat2Entrez),
                                    function(goid) {
                                        setdiff(cat2Entrez[[goid]],
                                                curCat2Entrez[[goid]])
                                    })

        ## White balls removed from urn by conditioning
        numSelectedRemoved <- sapply(cat2RemovedEntrez,
                                     function(x) sum(selected %in% x))

        ## Black balls removed from urn by conditioning
        numOtherRemoved <- listLen(cat2RemovedEntrez) - numSelectedRemoved
    } else {
        numSelectedRemoved <- rep(0, length(curCat2Entrez))
        numOtherRemoved <- numSelectedRemoved
    }

    ## Num white drawn (n11)
    numWdrawn <- sapply(curCat2Entrez, 
                        function(x) sum(selected %in% x))

    ## Num white
    numW <- listLen(curCat2Entrez)
    
    ## Num black
    numB <- (length(universeGeneIds(p)) - numOtherRemoved - numSelectedRemoved
             - numW)

    ## Num drawn
    numDrawn <- length(selected) - numSelectedRemoved

    n21 <- numW - numWdrawn
    n12 <- numDrawn - numWdrawn
    n22 <- numB - n12

    odds_ratio <-  (numWdrawn * n22) / (n12 * n21)

    expected <- (numWdrawn + n12) * (numWdrawn + n21)
    expected <- expected / (numWdrawn + n21 + n21 + n22)

    if (testDirection(p) == "over") {
        ## take the -1 because we want evidence for as extreme or more
        pvals <- phyper(numWdrawn - 1, numW, numB,
                        numDrawn, lower.tail=FALSE)
    } else {
        pvals <- phyper(numWdrawn, numW, numB,
                        numDrawn, lower.tail=TRUE)
    }
    list(p=pvals, odds=odds_ratio, expected=expected)
}
