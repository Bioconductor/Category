setClass("DatPkg",
         contains="VIRTUAL",
         representation=representation(
           name="character"))

setClass("AffyDatPkg", contains="DatPkg")
setClass("YeastDatPkg", contains="DatPkg")
setClass("ArabidopsisDatPkg", contains="DatPkg")

## For hummanLLMapping and similar
## setClass("OrganismMappingDatPkg", contains="DatPkg")

setClass("Org.XX.egDatPkg", contains="DatPkg")

DatPkgFactory <- function(chip) {
    
    strMatch <- function(pat, s) length(grep(pat, s)) > 0

    if(strMatch(".db$",chip)) chip<- sub(".db","",chip)
    
    if (chip == "UNKNOWN")
      return(new("AffyDatPkg", name=chip))

    pkg = paste(chip,".db",sep="")
    
    ##Use standardized schema names to decide
    if(require(pkg, character.only = TRUE)){
        conn <- do.call(paste(chip, "_dbconn", sep=""), list())
        schema <- dbmeta(conn, "DBSCHEMA")
        
        if(schema == "YEAST_DB" || schema == "YEASTCHIP_DB")
          pkg <- new("YeastDatPkg", name=chip)
        else if( schema == "ARABIDOPSIS_DB" || schema == "ARABIDOPSISCHIP_DB" )
          pkg <- new("ArabidopsisDatPkg", name=chip)
        else if( strMatch("CHIP_DB$", schema)){
            pkg <- new("AffyDatPkg", name=chip)}
        else { ##Otherwise its an ordinary org package
            pkg <- new("Org.XX.egDatPkg", name=chip)
        }
        return(pkg)
    }else stop(paste("Required annotation package", chip, "is not available.",sep=" "))
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

setClass("KEGGHyperGResult",
         contains="HyperGResult")

setClass("PFAMHyperGResult",
         contains="HyperGResult")

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




### DS: Similar structures for 'linear model'-based tests. minSize is
### needed mostly for conditional tests, becausewe don't want to leave
### out sub-genesets that are too small, even if they are significant
### (or do we?)


setClass("LinearMParams",
         contains="VIRTUAL",
         representation =
         representation(geneStats="numeric",
                        universeGeneIds="ANY",
                        annotation="character",
                        datPkg="DatPkg",
                        cateogrySubsetIds="ANY",
                        categoryName="character",
                        pvalueCutoff="numeric",
                        minSize="integer",
                        testDirection="character"), ## less, greater, two-sided?
         prototype=
         prototype(pvalueCutoff=0.01,
                   testDirection="up",
                   minSize=5L,
                   datPkg=DatPkgFactory("UNKNOWN")))

setClass("ChrMapLinearMParams",
         contains="LinearMParams",
         representation =
         representation(chrGraph="graph",
                        conditional="logical"),
         prototype =
         prototype(categoryName="ChrMap",
                   chrGraph=new("graphNEL", edgemode="directed")))





setClass("LinearMResultBase",
         representation(annotation="character",
                        geneIds="ANY",
                        testName="character",
                        pvalueCutoff="numeric",
                        minSize="integer",
                        testDirection="character"),
         contains="VIRTUAL",
         prototype=prototype(pvalueCutoff=0.01))


setClass("LinearMResult",
         contains="LinearMResultBase",
         representation=
         representation(pvalues="numeric",
                        effectSize="numeric",
                        catToGeneId="list"))


setClass("ChrMapLinearMResult",
         contains="LinearMResult", ## FIXME: is hyperG version correct?
         representation=
         representation(pvalue.order="integer",
                        conditional="logical",
                        chrGraph="graph"),
         prototype=
         prototype(chrGraph=new("graphNEL", edgemode="directed")))

