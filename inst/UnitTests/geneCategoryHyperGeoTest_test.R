makeSimpleGeneGoHyperGeoTestParams <- function() {
    set.seed(344)
    probeIds <- ls(hgu95av2LOCUSID)
    randProbeIds <- sample(probeIds, 500)
##     entrezUniverse <- unlist(mget(randProbeIds, hgu95av2LOCUSID,
##                                   ifnotfound=NA))
    ## This is "wrong", should unlist, but the code
    ## should catch/correct it.  The right way is above.
    entrezUniverse <- mget(randProbeIds, hgu95av2LOCUSID,
                           ifnotfound=NA)

    entrezUniverse <- entrezUniverse[!is.na(entrezUniverse)]
    selectedEntrezIds <- sample(entrezUniverse, 30)
    params <- new("GeneGoHyperGeoTestParams",
                  geneIds=selectedEntrezIds, 
                  universeGeneIds=entrezUniverse,
                  annotation="hgu95av2", 
                  ontology="BP",
                  pvalue.cutoff=0.05,
                  conditional=FALSE,
                  test.direction="over")
    params
}
    

test_basic_regression <- function() {
    p <- makeSimpleGeneGoHyperGeoTestParams()
    res <- geneCategoryHyperGeoTest(p)
    checkEquals(18, sum(pvalues(res) < res@pvalue.cutoff))

    pvals <- round(c(0.01596240, 0.01825930, 0.02194761), 3)
    names(pvals) <- c("GO:0043170", "GO:0044265", "GO:0009057") 
    checkEquals(pvals, round(pvalues(res)[1:3], 3))
    
    gcounts <- c(13, 4, 4)
    names(gcounts) <- c("GO:0043170", "GO:0044265", "GO:0009057")
    checkEquals(gcounts, geneCounts(res)[1:3])

    ucounts <- c(134, 19, 20)
    names(ucounts) <- c("GO:0043170", "GO:0044265", "GO:0009057")
    checkEquals(ucounts, universeCounts(res)[1:3])
    
    checkEquals(381, universeMappedCount(res))
    checkEquals(22, geneMappedCount(res))
    checkEquals("hgu95av2", annotation(res))
    checkEquals(c("GO", "BP"), testName(res))
}
