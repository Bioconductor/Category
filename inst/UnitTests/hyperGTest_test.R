library("hgu95av2")
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

    pvals <- round(c(0.001, 0.001, 0.002), 3)
    names(pvals) <- c("GO:0051247", "GO:0009891", "GO:0044267") 
    checkEquals(pvals, round(pvalues(res)[3:5], 3))
    
    gcounts <- c(3, 3, 13)
    names(gcounts) <- names(pvals)
    checkEquals(gcounts, geneCounts(res)[3:5])

    ucounts <- c(4, 4, 95)
    names(ucounts) <- names(pvals)
    checkEquals(ucounts, universeCounts(res)[3:5])
    
    checkEquals(390, universeMappedCount(res))
    checkEquals(25, geneMappedCount(res))
    checkEquals("hgu95av2", annotation(res))
    checkEquals(c("GO", "BP"), testName(res))
}


test_basic_regression_YEAST <- function() {
    p <- makeSimpleGOHyperGParams(chip="YEAST")
    res <- hyperGTest(p)
    checkEquals(2, sum(pvalues(res) < pvalueCutoff(res)))

    pvals <- round(c(0.06, 0.06, 0.06), 3)
    names(pvals) <- c("GO:0045332", "GO:0007009", "GO:0045490")
    checkEquals(pvals, round(pvalues(res)[3:5], 3))
    
    gcounts <- c(1, 1, 1)
    names(gcounts) <- names(pvals)
    checkEquals(gcounts, geneCounts(res)[3:5])

    ucounts <- c(1, 1, 1)
    names(ucounts) <- names(pvals)
    checkEquals(ucounts, universeCounts(res)[3:5])
    
    checkEquals(500, universeMappedCount(res))
    checkEquals(30, geneMappedCount(res))
    checkEquals("YEAST", annotation(res))
    checkEquals(c("GO", "BP"), testName(res))
}
