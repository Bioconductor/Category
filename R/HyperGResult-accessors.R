## Accessor methods for HyperGResult class

setMethod("annotation", signature(object="HyperGResultBase"),
          function(object) object@annotation)

setMethod("pvalues", signature(r="HyperGResult"),
          function(r) r@pvalues)

setMethod("geneCounts", signature(r="HyperGResult"),
          function(r) r@geneCounts)

setMethod("universeCounts", signature(r="HyperGResult"),
          function(r) r@universeCounts)

setMethod("universeMappedCount", signature(r="HyperGResult"),
           function(r) length(unique(unlist(r@catToGeneId))))

setMethod("geneMappedCount", signature(r="HyperGResultBase"),
           function(r) length(r@geneIds))

## generic "annotation" defined in Biobase

setMethod("geneIds", signature(r="HyperGResultBase"),
          function(r) r@geneIds)

setMethod("testName", signature(r="HyperGResultBase"),
          function(r) r@testName)

setMethod("pvalueCutoff", signature(r="HyperGResultBase"),
          function(r) r@pvalueCutoff)

setMethod("testDirection", signature(r="HyperGResultBase"),
          function(r) r@testDirection)

setMethod("description",
          signature(object="HyperGResultBase"),
          function(object) {
              desc <- paste("Gene to %s Category Test for %s Representation",
                            "Test Result")
              desc <- sprintf(desc, paste(testName(object), collapse=" "),
                              testDirection(object))
              desc
          })

setMethod("isConditional", "HyperGResultBase",
          function(r) FALSE)
          
