setMethod("geneCategoryHyperGeoTest",
          signature(p="GeneCategoryHyperGeoTestParams"), 
          function(p) {
              origGeneIds <- p@geneIds
              p@universeGeneIds <- universeBuilder(p)
              selected <- intersect(p@geneIds, p@universeGeneIds)
              p@geneIds <- selected
              cat2Entrez <- categoryToEntrezBuilder(p)
              numFound <- sapply(cat2Entrez, function(x) sum(selected %in% x))
              numDrawn <- length(selected)
              ## num white in urn
              numAtCat <- sapply(cat2Entrez, length)
              ## num black in urn
              numNotAtCat <- length(p@universeGeneIds) - numAtCat
              ## take the -1 because we want evidence for as extreme or more.
              pvals <- phyper(numFound-1, numAtCat, numNotAtCat, numDrawn,
                              lower.tail=FALSE)
              ord <- order(pvals)
              new("GeneCategoryHyperGeoTestResult",
                  pvalues=pvals[ord],
                  geneCounts=numFound[ord],
                  universeCounts=numAtCat[ord],
                  catToGeneId=cat2Entrez,
                  annotation=p@annotation,
                  geneIds=p@geneIds,
                  testName=categoryName(p),
                  pvalue.cutoff=p@pvalue.cutoff)
          })


geneGoHyperGeoTest <- function(entrezGeneIds, lib, ontology, universe=NULL)
{
    if (missing(universe) || is.null(universe))
      universe <- character(0)
    params <- new("GeneGoHyperGeoTestParams",
                  geneIds=entrezGeneIds,
                  universeGeneIds=universe,
                  annotation=lib,
                  ontology=ontology)
    geneCategoryHyperGeoTest(params)
}

resultToGOHyperG <- function(res, origIds) {
    go2Affy <- mget(names(pvalues(res)),
                    getDataEnv("GO2ALLPROBES", annotation(res)))
    list(pvalues=pvalues(res),
         goCounts=universeCounts(res),
         intCounts=geneCounts(res),
         numLL=universeMappedCount(res),
         numInt=geneMappedCount(res),
         chip=annotation(res),
         intLLs=origIds,
         go2Affy=go2Affy)
}

GOHyperG <- function(x, lib, what="MF", universe=NULL)
{
    if (missing(universe) || is.null(universe))
      universe <- character(0)
    params <- new("GeneGoHyperGeoTestParams",
                  geneIds=x,
                  universeGeneIds=universe,
                  annotation=lib,
                  ontology=what)
    if (missing(lib) || !is.character(lib))
      stop("argument ", sQuote("lib"), " must be character")

    res <- geneCategoryHyperGeoTest(params)
    resultToGOHyperG(res, origIds=x)
}


geneKeggHyperGeoTest <- function(entrezGeneIds, lib, universe=NULL)
{
    if (missing(universe) || is.null(universe))
      universe <- character(0)
    params <- new("GeneKeggHyperGeoTestParams",
                  geneIds=entrezGeneIds,
                  universeGeneIds=universe,
                  annotation=lib)
    geneCategoryHyperGeoTest(params)
}
