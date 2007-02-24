library("Category")
library("YEAST")

test_KEGG1 <- function() {
    set.seed(434)    
    allYeast <- ls(YEASTCHR)
    selGenes <- sample(allYeast, 80)
    kp <- new("KEGGHyperGParams",
              geneIds=selGenes,
              annotation="YEAST")
    ans <- hyperGTest(kp)
    checkEquals("KEGG", testName(ans))
    checkEquals(14, length(geneIds(ans)))
    checkEquals(27, length(universeCounts(ans)))
    checkEquals(27, length(geneCounts(ans)))
}
