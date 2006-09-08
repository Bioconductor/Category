setClass("GeneCategoryHyperGeoTestParams", 
         representation(geneIds="ANY",
                        universeGeneIds="ANY",
                        annotation="character",
                        cateogrySubsetIds="ANY",
                        categoryName="character",
                        pvalueCutoff="numeric",
                        testDirection="character"),
         prototype=prototype(
           pvalueCutoff=0.01,
           testDirection="over"
           ),  ## FIXME: add validity check
         contains="VIRTUAL")


setClass("GeneGoHyperGeoTestParams",
         representation(ontology="character",
                        conditional="logical"),
         contains="GeneCategoryHyperGeoTestParams",
         prototype=prototype(categoryName="GO",
           conditional=FALSE))


setClass("GeneKeggHyperGeoTestParams",
         contains="GeneCategoryHyperGeoTestParams",
         prototype=prototype(categoryName="KEGG"))


setClass("GeneCategoryHyperGeoTestResultBase",
         representation(annotation="character",
                        geneIds="ANY",
                        testName="character",
                        pvalueCutoff="numeric",
                        testDirection="character"),
         contains="VIRTUAL",
         prototype=prototype(pvalueCutoff=0.01))


setClass("GeneCategoryHyperGeoTestResult",
         contains="GeneCategoryHyperGeoTestResultBase",
         representation=representation(pvalues="numeric",
           geneCounts="integer",
           universeCounts="integer",
           catToGeneId="list"))

