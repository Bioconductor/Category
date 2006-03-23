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
                  universeMappedCount=length(p@universeGeneIds),
                  geneMappedCount=length(p@geneIds),
                  annotation=p@annotation,
                  geneIds=origGeneIds,
                  testName=p@categoryName)
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


GOHyperG <- function(x, lib="hgu95av2", what="MF", universe=NULL)
{
    if (missing(universe) || is.null(universe))
      universe <- character(0)
    params <- new("GeneGoHyperGeoTestParams",
                  geneIds=x,
                  universeGeneIds=universe,
                  annotation=lib,
                  ontology=what)

    res <- geneCategoryHyperGeoTest(params)
    go2Affy <- mget(names(res@pvalues), getDataEnv("GO2ALLPROBES", lib))
    list(pvalues=res@pvalues,
         goCounts=res@universeCounts,
         intCounts=res@geneCounts,
         numLL=res@universeMappedCount,
         numInt=res@geneMappedCount,
         chip=res@annotation,
         intLLs=res@geneIds,
         go2Affy=go2Affy)
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
