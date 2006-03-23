setClass("GeneCategoryHyperGeoTestParams", 
         representation(geneIds="ANY",
                        universeGeneIds="ANY",
                        annotation="character",
                        cateogrySubsetIds="ANY",
                        categoryName="character"),
         contains="VIRTUAL")


setClass("GeneGoHyperGeoTestParams", representation(ontology="character"),
         contains="GeneCategoryHyperGeoTestParams",
         prototype=prototype(categoryName="GO"))


setClass("GeneKeggHyperGeoTestParams", contains="GeneCategoryHyperGeoTestParams",
         prototype=prototype(categoryName="KEGG"))


setClass("GeneCategoryHyperGeoTestResult",
         representation(pvalues="numeric",
                        geneCounts="integer",
                        universeCounts="integer",
                        universeMappedCount="integer",
                        geneMappedCount="integer",
                        annotation="character",
                        geneIds="ANY",
                        testName="character"))
