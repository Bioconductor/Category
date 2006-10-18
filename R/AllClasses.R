setClass("DatPkg",
         contains="VIRTUAL",
         representation=representation(name="character"))

setClass("AffyDatPkg", contains="DatPkg")
setClass("YeastDatPkg", contains="DatPkg")
## For hummanLLMapping and similar
setClass("OrganismMappingDatPkg", contains="DatPkg")

DatPkgFactory <- function(pkgName) {
    if (length(grep("YEAST", pkgName)) > 0)
      pkg <- new("YeastDatPkg", name=pkgName)
    else if (length(grep("Mapping", pkgName)) > 0)
      pkg <- new("OrganismMappingDatPkg", name=pkgName)
    else
      pkg <- new("AffyDatPkg", name=pkgName)
    pkg
}


setClass("HyperGParams", 
         representation(geneIds="ANY",
                        universeGeneIds="ANY",
                        annotation="character",
                        datPkg="DatPkg",
                        cateogrySubsetIds="ANY",
                        categoryName="character",
                        pvalueCutoff="numeric",
                        testDirection="character"),
         prototype=prototype(
           pvalueCutoff=0.01,
           testDirection="over",
           datPkg=DatPkgFactory("UNKNOWN")
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

setClass("PFAMHyperGParams",
         contains="HyperGParams",
         prototype=prototype(categoryName="PFAM"))


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


