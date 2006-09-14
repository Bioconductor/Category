setClass("HyperGParams", 
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


setClass("GOHyperGParams",
         representation(ontology="character",
                        conditional="logical"),
         contains="HyperGParams",
         prototype=prototype(categoryName="GO",
           conditional=FALSE))


setClass("KEGGHyperGParams",
         contains="HyperGParams",
         prototype=prototype(categoryName="KEGG"))


setClass("HyperGResultBase",
         representation(annotation="character",
                         geneIds="ANY",
                        testName="character",
                        pvalueCutoff="numeric",
                        testDirection="character"),
         contains="VIRTUAL",
         prototype=prototype(pvalueCutoff=0.01))


setClass("HyperGResult",
         contains="HyperGResultBase",
         representation=representation(
           pvalues="numeric",
           oddsRatios="numeric",
           expectedCounts="numeric",
           geneCounts="integer",
           universeCounts="integer",
           catToGeneId="list"))

