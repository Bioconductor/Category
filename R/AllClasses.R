setClass("DatPkg",
         contains="VIRTUAL",
         representation=representation(name="character"))

setClass("DBPkg",
         contains=c("DatPkg", "VIRTUAL",
         representation=representation(getdb="function"))

setClass("AffyDatPkg", contains="DatPkg")
setClass("AffyDBPkg", contains="DBPkg")

setClass("YeastDatPkg", contains="DatPkg")
setClass("YeastDBPkg", contains="DBPkg")

## For hummanLLMapping and similar
setClass("OrganismMappingDatPkg", contains="DatPkg")
setClass("OrganismMappingDBPkg", contains="DBPkg")


DatPkgFactory <- function(pkgName) {

    strMatch <- function(pat, s) length(grep(p, s)) > 0
    isDbPkg <- function(p) length(grep("db$", p)) > 0

    havePkg <- require(pkgName, character.only=TRUE)
    if (!havePkg)
      stop("the ", pkgName, " package was not found.")

    ## XXX: ugly name-based computations ahead
    if (!isDbPkg(pkgName)) {
        if (strMatch("YEAST", pkgName))
          pkg <- new("YeastDatPkg", name=pkgName)
        else if (strMatch("Mapping", pkgName))
          pkg <- new("OrganismMappingDatPkg", name=pkgName)
        else
          pkg <- new("AffyDatPkg", name=pkgName)
    } else { ## we have a db-based package
        pkgns <- getNamespace(pkgName)
        dbgetter <- get("getDb", envir=pkgns)
        if (strMatch("YEAST", pkgName))
          pkg <- new("YeastDBPkg", name=pkgName, getdb=dbgetter)
        else if (strMatch("Mapping", pkgName))
          pkg <- new("OrganismMappingDBPkg", name=pkgName,
                     getdb=dbgetter)
        else
          pkg <- new("AffyDBPkg", name=pkgName, getdb=dbgetter)
    }
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


