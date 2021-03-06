### Accessor methods for HyperGResultBase class

setMethod("annotation", signature(object="HyperGResultBase"),
          function(object) object@annotation)

setMethod("geneCounts", signature(r="HyperGResultBase"),
          function(r) {
              sapply(geneIdUniverse(r), function(x) {
                  sum(geneIds(r) %in% x)
              })
          })

setMethod("universeCounts", signature(r="HyperGResultBase"),
          function(r) {
              univ <- geneIdUniverse(r)
              ans <- listLen(univ)
              names(ans) <- names(univ)
              ans
          })

setMethod("universeMappedCount", signature(r="HyperGResultBase"),
           function(r) length(unique(unlist(geneIdUniverse(r)))))

setMethod("geneMappedCount", signature(r="HyperGResultBase"),
           function(r) length(geneIds(r)))

setMethod("geneIds", signature(object="HyperGResultBase"),
          function(object, ...) object@geneIds)

setMethod("geneIdsByCategory", signature(r="HyperGResultBase"),
          function(r, catids=NULL) {
              ans <- geneIdUniverse(r)
              if (!missing(catids) && !is.null(catids))
                ans <- ans[catids]
              lapply(ans, intersect, geneIds(r))
          })

setMethod("sigCategories", signature(r="HyperGResultBase"),
          function(r, p) {
              if (missing(p))
                p <- pvalueCutoff(r)
              pv <- pvalues(r)
              names(pv[pv < p])
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
              cond <- "Conditional"
              if (!conditional(object))
                cond <- ""
              desc <- paste("Gene to %s", cond, "test for %s-representation")
              desc <- sprintf(desc, paste(testName(object), collapse=" "),
                              testDirection(object))
              desc
          })

setMethod("conditional", "HyperGResultBase",
          function(r) {
              if (!("conditional" %in% slotNames(r)))
                FALSE
              else
                r@conditional
          })



### Accessor methods for HyperGResult class

setMethod("pvalues", signature(r="HyperGResult"),
          function(r) r@pvalues)

setMethod("oddsRatios", signature(r="HyperGResult"),
          function(r) r@oddsRatios)

setMethod("expectedCounts", signature(r="HyperGResult"),
          function(r) r@expectedCounts)

## generic "annotation" defined in Biobase

setMethod("geneIdUniverse", signature(r="HyperGResult"),
          function(r, cond=TRUE) r@catToGeneId)
