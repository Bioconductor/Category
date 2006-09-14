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
              numFound <- sapply(cat2Entrez, function(x) sum(selected %in% x))
              numDrawn <- length(selected)
              ## num white in urn
              numAtCat <- sapply(cat2Entrez, length)
              ## num black in urn
              numNotAtCat <- length(universeGeneIds(p)) - numAtCat
              ## take the -1 because we want evidence for as extreme or more.
              pvals <- phyper(numFound-1, numAtCat, numNotAtCat, numDrawn,
                              lower.tail=FALSE)
              ord <- order(pvals)
              new("HyperGResult",
                  pvalues=pvals[ord],
                  geneCounts=numFound[ord],
                  universeCounts=numAtCat[ord],
                  catToGeneId=cat2Entrez,
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
