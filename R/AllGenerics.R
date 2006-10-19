setGeneric("hyperGTest", 
           function(p) standardGeneric("hyperGTest"),
           valueClass="HyperGResultBase")

setGeneric("categoryToEntrezBuilder", 
           function(p) standardGeneric("categoryToEntrezBuilder"))

setGeneric("universeBuilder", 
           function(p) standardGeneric("universeBuilder"))


## Accessors for HyperGResult objects
setGeneric("pvalues", function(r) standardGeneric("pvalues"))

setGeneric("geneCounts", function(r) standardGeneric("geneCounts"))

setGeneric("universeCounts", function(r) standardGeneric("universeCounts"))

setGeneric("universeMappedCount",
           function(r) standardGeneric("universeMappedCount"))

setGeneric("geneMappedCount",
           function(r) standardGeneric("geneMappedCount"))

## generic "annotation" defined in Biobase

setGeneric("geneIds", function(r) standardGeneric("geneIds"))
setGeneric("geneIds<-", function(r, value) standardGeneric("geneIds<-"))

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

setGeneric("ontology", function(r) standardGeneric("ontology"))

setGeneric("conditional", function(r) standardGeneric("conditional"))

setGeneric("categoryName<-", function(r, value) standardGeneric("categoryName<-"))

setGeneric("universeGeneIds<-",
           function(r, value) standardGeneric("universeGeneIds<-"))

setGeneric("ontology<-", function(r, value) standardGeneric("ontology<-"))

setGeneric("conditional<-", function(r, value) standardGeneric("conditional<-"))

setGeneric("categorySubsetIds", function(r) standardGeneric("categorySubsetIds"))

setGeneric("categorySubsetIds<-", function(r, value) standardGeneric("categorySubsetIds<-"))

setGeneric("isConditional", function(r) standardGeneric("isConditional"))


## accessors for DatPkg
setGeneric("ID2GO", function(p) standardGeneric("ID2GO"))
setGeneric("ID2EntrezID", function(p) standardGeneric("ID2EntrezID"))
setGeneric("GO2AllProbes", signature=c("p"),
           function(p, ontology) standardGeneric("GO2AllProbes"))
