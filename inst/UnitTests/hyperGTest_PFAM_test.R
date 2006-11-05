library("Category")
library("YEAST")

test_PFAM1 <- function() {
    set.seed(434)    
    allYeast <- ls(YEASTCHR)
    selGenes <- sample(allYeast, 80)
    pp <- new("PFAMHyperGParams",
              geneIds=selGenes,
              annotation="YEAST")
    ans <- hyperGTest(pp)
    checkEquals("PFAM", testName(ans))
    checkEquals(38, length(geneIds(ans)))
    checkEquals(50, length(universeCounts(ans)))
}
