library("hgu95av2.db")
library("YEAST")

makeSimpleGOHyperGParams <- function(chip=c("hgu95av2", "YEAST")) {
    set.seed(344)
    chip <- match.arg(chip)
    if (chip == "hgu95av2") {
        probeIds <- ls(hgu95av2ENTREZID)
    } else if (chip == "YEAST") {
        probeIds <- ls(YEASTCHR)
    }
    randProbeIds <- sample(probeIds, 500)
    if (chip == "YEAST") {
        entrezUniverse <- randProbeIds
    } else {
        entrezUniverse <- mget(randProbeIds, hgu95av2ENTREZID,
                               ifnotfound=NA)
    }
    entrezUniverse <- entrezUniverse[!is.na(entrezUniverse)]
    selectedEntrezIds <- sample(entrezUniverse, 30)
    params <- new("GOHyperGParams",
                  geneIds=selectedEntrezIds, 
                  universeGeneIds=entrezUniverse,
                  annotation=chip, 
                  ontology="BP",
                  pvalueCutoff=0.05,
                  conditional=FALSE,
                  testDirection="over")
    params
}
    

test_basic_regression_hgu95av2 <- function() {
    p <- makeSimpleGOHyperGParams(chip="hgu95av2")
    res <- hyperGTest(p)
    checkEquals(22, sum(pvalues(res) < pvalueCutoff(res)))

    ## This is a fragile test.  It depends on the annotation
    ## data.  Choose pvalues to compare carefully since many
    ## terms will have same pvalue and you don't want to fail
    ## just because the ordering is different.
    pvals <- round(c(0.126, 0.129, 0.138), 3)
    names(pvals) <- c("GO:0048522", "GO:0006968", "GO:0048518")
    checkEquals(pvals, round(pvalues(res)[87:89], 3))
    
    gcounts <- c(4, 2, 4)
    names(gcounts) <- names(pvals)
    checkEquals(gcounts, geneCounts(res)[87:89])

    ucounts <- c(31, 10, 32)
    names(ucounts) <- names(pvals)
    checkEquals(ucounts, universeCounts(res)[87:89])
    
    checkEquals(390, universeMappedCount(res))
    checkEquals(25, geneMappedCount(res))
    checkEquals("hgu95av2", annotation(res))
    checkEquals(c("GO", "BP"), testName(res))
}


test_basic_regression_YEAST <- function() {
    p <- makeSimpleGOHyperGParams(chip="YEAST")
    res <- hyperGTest(p)
    checkEquals(2, sum(pvalues(res) < pvalueCutoff(res)))

    pvals <- round(c(0.927, 0.929, 0.956 ), 3)
    names(pvals) <- c("GO:0044237", "GO:0006412", "GO:0009059")
    checkEquals(pvals, round(pvalues(res)[201:203], 3))
    
    gcounts <- c(12, 1, 1)
    names(gcounts) <- names(pvals)
    checkEquals(gcounts, geneCounts(res)[201:203])

    ucounts <- c(256, 41, 48)
    names(ucounts) <- names(pvals)
    checkEquals(ucounts, universeCounts(res)[201:203])
    
    checkEquals(500, universeMappedCount(res))
    checkEquals(30, geneMappedCount(res))
    checkEquals("YEAST", annotation(res))
    checkEquals(c("GO", "BP"), testName(res))
}
