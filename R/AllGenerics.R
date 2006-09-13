setGeneric("geneCategoryHyperGeoTest", 
           function(p) standardGeneric("geneCategoryHyperGeoTest"),
           valueClass="HyperGResult")

setGeneric("categoryToEntrezBuilder", 
           function(p) standardGeneric("categoryToEntrezBuilder"))

setGeneric("universeBuilder", 
           function(p) standardGeneric("universeBuilder"))

setGeneric("categoryName", function(p) standardGeneric("categoryName"))


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

setGeneric("testName", function(r) standardGeneric("testName"))

setGeneric("pvalueCutoff", function(r) standardGeneric("pvalueCutoff"))

setGeneric("testDirection", function(r) standardGeneric("testDirection"))
