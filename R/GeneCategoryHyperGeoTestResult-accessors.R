## Accessor methods for GeneCategoryHyperGeoTestResult class

setMethod("annotation", signature(object="GeneCategoryHyperGeoTestResultBase"),
          function(object) object@annotation)

setMethod("pvalues", signature(r="GeneCategoryHyperGeoTestResult"),
          function(r) r@pvalues)

setMethod("geneCounts", signature(r="GeneCategoryHyperGeoTestResult"),
          function(r) r@geneCounts)

setMethod("universeCounts", signature(r="GeneCategoryHyperGeoTestResult"),
          function(r) r@universeCounts)

setMethod("universeMappedCount", signature(r="GeneCategoryHyperGeoTestResult"),
           function(r) length(unique(unlist(r@catToGeneId))))

setMethod("geneMappedCount", signature(r="GeneCategoryHyperGeoTestResultBase"),
           function(r) length(r@geneIds))

## generic "annotation" defined in Biobase

setMethod("geneIds", signature(r="GeneCategoryHyperGeoTestResultBase"),
          function(r) r@geneIds)

setMethod("testName", signature(r="GeneCategoryHyperGeoTestResultBase"),
          function(r) r@testName)

setMethod("pvalueCutoff", signature(r="GeneCategoryHyperGeoTestResultBase"),
          function(r) r@pvalue.cutoff)

