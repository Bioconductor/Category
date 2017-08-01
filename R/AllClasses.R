setClass("DatPkg",
         contains="VIRTUAL",
         representation=representation(
           name="character"))

##Annotion package DatPkgs
setClass("AffyDatPkg", contains="DatPkg")
setClass("YeastDatPkg", contains="DatPkg")
setClass("ArabidopsisDatPkg", contains="DatPkg")
setClass("Org.XX.egDatPkg", contains="DatPkg")

## Other types of DatPkgs
setClass("GeneSetCollectionDatPkg", contains="DatPkg",
         representation=representation(
           geneSetCollection="GeneSetCollection"))

## OBO gene set collection DatPkgs
setClass("OBOCollectionDatPkg", contains="DatPkg",
         representation=representation(
           oboCollection="OBOCollection",
           oboGraph="graph",
           geneSetCollection="GeneSetCollection"))

## These generics needed for AllClasses
setGeneric("configureDatPkg",
           function(annotation, ...) standardGeneric("configureDatPkg"))

setGeneric("DatPkgFactory",
           function(chip) standardGeneric("DatPkgFactory"),
           useAsDefault=function(chip) new("AffyDatPkg", name="UNKNOWN"))

setClass("HyperGParams",
         contains="VIRTUAL",
         representation=representation(geneIds="ANY",
           universeGeneIds="ANY",
           annotation="character",
           datPkg="DatPkg",
           categorySubsetIds="ANY",
           categoryName="character",
           pvalueCutoff="numeric",
           testDirection="character"),
         prototype=prototype(
           pvalueCutoff=0.01,
           testDirection="over",
           datPkg=DatPkgFactory()
           ))


setClass("GOHyperGParams",
         representation(ontology="character",
                        conditional="logical",
                        orCutoff="numeric",
                        minSizeCutoff="numeric",
                        maxSizeCutoff="numeric"),
         contains="HyperGParams",
         prototype=prototype(
           categoryName="GO",
           conditional=FALSE,
           orCutoff=1,
           minSizeCutoff=0,
           maxSizeCutoff=Inf,
           annotation="GO"))


setClass("KEGGHyperGParams",
         contains="HyperGParams",
         prototype=prototype(
           categoryName="KEGG",
           annotation="KEGG"))

setClass("PFAMHyperGParams",
         contains="HyperGParams",
         prototype=prototype(categoryName="PFAM"))

## this probably should not go here
setClass("GeneSetCollectionAnnotation", contains="character")

.GeneSetCollectionAnnotation <- function(annotation)
      new("GeneSetCollectionAnnotation", annotation)

## class for OBO gene set testing
setClass("OBOHyperGParams",
         contains="HyperGParams",
         representation(conditional="logical",
                        orCutoff="numeric",
                        minSizeCutoff="numeric",
                        maxSizeCutoff="numeric"),
         prototype=prototype(
           categoryName="OBO",
           conditional=FALSE,
           orCutoff=1,
           minSizeCutoff=0,
           maxSizeCutoff=Inf,
           annotation=.GeneSetCollectionAnnotation("OBO")))

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
           catToGeneId="list",
           organism="character"))

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

## this is for OBO
setClass("OBOHyperGResult",
         contains="HyperGResultBase",
         representation=representation(
           goDag="graph",
           gscDescriptions="character",
           pvalue.order="integer",
           conditional="logical"),
         prototype=prototype(
           testName="GO",
           pvalueCutoff=0.01,
           goDag=new("graphNEL")))



### DS: Similar structures for 'linear model'-based tests. minSize is
### needed mostly for conditional tests, because we don't want to leave
### out sub-genesets that are too small, even if they are significant
### (or do we?)

### ML: Pushed the ChrMap stuff up to the general classes; not sure
### why we even need a class for each type of category? Seems too
### formal for something that can be so arbitrary. Every type of
### category (GO, KEGG, chromosome bands, ...) fits into a graph.

### ML: Current thoughts are that every Params object should have a
### GeneSetCollection that refers (through its CollectionType) to an
### ontology. From that, a graph can be automatically derived. This
### favors composition over inheritance.

setClass("LinearMParams",
         representation =
         representation(geneStats="numeric", # this needs to be named
                        universeGeneIds="ANY", # not used, probably should be
                        annotation="character", # just carried through...
                        datPkg="DatPkg", # not currently used
                        categorySubsetIds="ANY", # not used, probably should be
                        categoryName="character",
                        pvalueCutoff="numeric",
                        minSize="integer",
                        testDirection="character",## less, greater, two-sided?
#### FIXME: could the graph be derived from the 'CollectionType's in the
#### 'GeneSetCollection', using info from e.g. GO.db?
                        graph = "graph",
                        conditional="logical",
                        ## instead of putting attributes on graph
                        gsc = "GeneSetCollection"), 
         prototype=
         prototype(pvalueCutoff=0.01,
                   testDirection="up",
                   minSize=5L,
                   datPkg=DatPkgFactory(),
                   conditional = FALSE,
                   graph = new("graphNEL", edgemode = "directed")),
         validity = function(object) {
           if (length(object@geneStats) && is.null(names(object@geneStats)))
             "'geneStats' must have names that match those in the gene sets"
           else NULL
         })

setClass("ChrMapLinearMParams",
         contains="LinearMParams",
         prototype =
         prototype(categoryName="ChrMap"))





setClass("LinearMResultBase",
         representation(annotation="character",
                        geneIds="ANY", # FIXME: should this be universeGeneIds?
                        testName="character",
                        pvalueCutoff="numeric",
                        minSize="integer",
                        testDirection="character",
                        conditional="logical",
                        graph = "graph",
                        gsc = "GeneSetCollection"),
         contains="VIRTUAL",
         prototype=prototype(pvalueCutoff=0.01))


setClass("LinearMResult",
         contains="LinearMResultBase",
         representation=
         representation(pvalues="numeric",
                        pvalue.order="integer",
                        effectSize="numeric"))


setClass("ChrMapLinearMResult",
         contains = "LinearMResult")
