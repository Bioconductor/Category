setClass("GeneCategoryHyperGeoTestParams", 
         representation(geneIds="ANY",
                        universeGeneIds="ANY",
                        annotation="character",
                        cateogrySubsetIds="ANY",
                        categoryName="character",
                        pvalue.cutoff="numeric",
                        test.direction="character"),
         prototype=prototype(
           pvalue.cutoff=0.01,
           test.direction="over"
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
                        pvalue.cutoff="numeric",
                        test.direction="character"),
         contains="VIRTUAL",
         prototype=prototype(pvalue.cutoff=0.01))


setClass("GeneCategoryHyperGeoTestResult",
         contains="GeneCategoryHyperGeoTestResultBase",
         representation=representation(pvalues="numeric",
           geneCounts="integer",
           universeCounts="integer",
           catToGeneId="list"))

