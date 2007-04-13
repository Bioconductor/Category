## Accessor methods for HyperGResult class

setMethod("annotation", signature(object="HyperGResultBase"),
          function(object) object@annotation)

setMethod("pvalues", signature(r="HyperGResult"),
          function(r) r@pvalues)

setMethod("oddsRatios", signature(r="HyperGResult"),
          function(r) r@oddsRatios)

setMethod("expectedCounts", signature(r="HyperGResult"),
          function(r) r@expectedCounts)

setMethod("geneCounts", signature(r="HyperGResult"),
          function(r) {
              sapply(r@catToGeneId, function(x) {
                  sum(geneIds(r) %in% x)
              })
          })

setMethod("universeCounts", signature(r="HyperGResult"),
          function(r) {
              ans <- listLen(r@catToGeneId)
              names(ans) <- names(r@catToGeneId)
              ans
          })

setMethod("universeMappedCount", signature(r="HyperGResult"),
           function(r) length(unique(unlist(r@catToGeneId))))

setMethod("geneMappedCount", signature(r="HyperGResultBase"),
           function(r) length(r@geneIds))

## generic "annotation" defined in Biobase

setMethod("geneIds", signature(r="HyperGResultBase"),
          function(r) r@geneIds)

setMethod("geneIdUniverse", signature(r="HyperGResult"),
          function(r) r@catToGeneId)

setMethod("geneIdsByCategory", signature(r="HyperGResultBase"),
          function(r, catids=NULL) {
              ans <- geneIdUniverse(r)
              if (!missing(catids) && !is.null(catids))
                ans <- ans[catids]
              lapply(ans, intersect, geneIds(r))
          })

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
          
