setGeneric("hyperGTest", 
           function(p) standardGeneric("hyperGTest"),
           valueClass="HyperGResultBase")

setGeneric("linearMTest", 
           function(p) standardGeneric("linearMTest"),
           valueClass="LinearMResultBase")

setGeneric("categoryToEntrezBuilder", 
           function(p) standardGeneric("categoryToEntrezBuilder"))

setGeneric("universeBuilder", 
           function(p) standardGeneric("universeBuilder"))


## Accessors for HyperGResult and LinearMResult objects
setGeneric("pvalues", function(r) standardGeneric("pvalues"))

setGeneric("effectSize", function(r) standardGeneric("effectSize"))

setGeneric("geneCounts", function(r) standardGeneric("geneCounts"))

setGeneric("universeCounts", function(r) standardGeneric("universeCounts"))

setGeneric("universeMappedCount",
           function(r) standardGeneric("universeMappedCount"))

setGeneric("geneMappedCount",
           function(r) standardGeneric("geneMappedCount"))

setGeneric("chrGraph",
           function(r) standardGeneric("chrGraph"))

setGeneric("geneIdUniverse", signature="r",
           function(r, cond=TRUE) standardGeneric("geneIdUniverse"))

setGeneric("condGeneIdUniverse",
           function(r) {
               .Defunct(msg=paste("use unconditional version with",
                             "'cond' argument"))
               standardGeneric("condGeneIdUniverse")
           })

## generic "annotation" defined in Biobase

#setGeneric("geneIds", function(r, ...) standardGeneric("geneIds"))
#setGeneric("geneIds<-", function(r, value) standardGeneric("geneIds<-"))

setGeneric("geneIdsByCategory", signature="r",
           function(r, catids = NULL) standardGeneric("geneIdsByCategory"))

setGeneric("sigCategories", signature="r",
           function(r, p) standardGeneric("sigCategories"))

setGeneric("testName", function(r) standardGeneric("testName"))

setGeneric("pvalueCutoff", function(r) standardGeneric("pvalueCutoff"))
setGeneric("pvalueCutoff<-", function(r, value) standardGeneric("pvalueCutoff<-"))

setGeneric("testDirection", function(r) standardGeneric("testDirection"))
setGeneric("testDirection<-",
           function(r, value) standardGeneric("testDirection<-"))

setGeneric("oddsRatios", function(r) standardGeneric("oddsRatios"))

setGeneric("expectedCounts",
           function(r) standardGeneric("expectedCounts"))


## accoessors for HyperGParams
setGeneric("categoryName", function(r) standardGeneric("categoryName"))

setGeneric("universeGeneIds", function(r) standardGeneric("universeGeneIds"))

#setGeneric("ontology", function(r) standardGeneric("ontology"))

setGeneric("conditional", function(r) standardGeneric("conditional"))

setGeneric("categoryName<-", function(r, value) standardGeneric("categoryName<-"))

setGeneric("universeGeneIds<-",
           function(r, value) standardGeneric("universeGeneIds<-"))

setGeneric("ontology<-", function(r, value) standardGeneric("ontology<-"))

setGeneric("conditional<-", function(r, value) standardGeneric("conditional<-"))

setGeneric("categorySubsetIds", function(r) standardGeneric("categorySubsetIds"))

setGeneric("categorySubsetIds<-", function(r, value) standardGeneric("categorySubsetIds<-"))

setGeneric("isConditional", function(r) {
    .Defunct("conditional")
    standardGeneric("isConditional")
    })

setGeneric("htmlReport", function(r, file="", append=FALSE, label="",
                                  digits=3, summary.args=NULL)
           standardGeneric("htmlReport"),
           signature=c("r"))

setGeneric("makeValidParams", function(object) {
           v <- standardGeneric("makeValidParams")
           if (class(object) != class(v))
             stop(paste("makeValidParams generic must return same class ",
                        "as its argument 'object'"))
           v
           })

## accessors for DatPkg
setGeneric("ID2GO", function(p) standardGeneric("ID2GO"))
setGeneric("ID2EntrezID", function(p) standardGeneric("ID2EntrezID"))
setGeneric("GO2AllProbes", signature=c("p"),
           function(p, ontology) standardGeneric("GO2AllProbes"))

## ChrBandTree accessors
setGeneric("lgeneIds", function(r, ...) standardGeneric("lgeneIds"))

setGeneric("childrenOf",
           function(g, n, ...) standardGeneric("childrenOf"))

setGeneric("parentOf",
           function(g, n, ...) standardGeneric("parentOf"))

setGeneric("allGeneIds",
           function(g, ...) standardGeneric("allGeneIds"))

setGeneric("treeLevels",
           function(g, ...) standardGeneric("treeLevels"))

setGeneric("level2nodes",
           function(g, level, ...) standardGeneric("level2nodes"))
