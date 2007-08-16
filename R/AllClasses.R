setClass("DatPkg",
         contains="VIRTUAL",
         representation=representation(
           name="character"))

setClass("AffyDatPkg", contains="DatPkg")
setClass("YeastDatPkg", contains="DatPkg")

## For hummanLLMapping and similar
setClass("OrganismMappingDatPkg", contains="DatPkg")

setClass("Org.XX.egDatPkg", contains="DatPkg")

DatPkgFactory <- function(chip) {

    strMatch <- function(pat, s) length(grep(pat, s)) > 0

    if (chip == "UNKNOWN")
      return(new("AffyDatPkg", name=chip))
    ## XXX: ugly name-based computations ahead
    if (strMatch("YEAST", chip))
      pkg <- new("YeastDatPkg", name=chip)
    else if (strMatch("LLMapping", chip))
      pkg <- new("OrganismMappingDatPkg", name=chip)
    else if (strMatch("^org\\.[a-zA-Z]+\\.eg\\.db$", chip))
      pkg <- new("Org.XX.egDatPkg", name=chip)
    else
      pkg <- new("AffyDatPkg", name=chip)
    pkg
}


setClass("HyperGParams",
         contains="VIRTUAL",
         representation=representation(geneIds="ANY",
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
           ))


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

setClass("ChrMapHyperGParams",
         contains="HyperGParams",
         representation=representation(
           chrGraph="graph",
           conditional="logical"),
         prototype=prototype(
           categoryName="ChrMap",
           chrGraph=new("graphNEL", edgemode="directed")))

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
           catToGeneId="list"))


setClass("ChrMapHyperGResult",
         contains="HyperGResultBase",
         representation=representation(
           pvalue.order="integer",
           conditional="logical",
           chrGraph="graph"),
         prototype=prototype(
           chrGraph=new("graphNEL", edgemode="directed")))

setClass("ChrBandTree",
         representation=representation(
           toParentGraph="graph",
           toChildGraph="graph",
           root="character",
           level2nodes="list"))

