setMethod("show", signature(object="HyperGResultBase"),
          function(object) {
              cat(description(object), "\n")
              nSig <- sum(pvalues(object) < object@pvalueCutoff[1])
              cat(length(pvalues(object)), testName(object), "ids tested ")
              cat("(", nSig, " have p < ", object@pvalueCutoff[1],
                  ")\n", sep="")
              cat("Selected gene set size:", geneMappedCount(object), "\n")
              cat("    Gene universe size:", universeMappedCount(object), "\n")
              cat("    Annotation package:", annotation(object), "\n")
          })


setMethod("show", signature(object="HyperGParams"),
          function(object) {
              cat("A", class(object), "instance\n")
              cat("  category:", object@categoryName, "\n")
              cat("annotation:", object@annotation, "\n")
          })


setMethod("show", signature(object="ChrBandTree"),
          function(object) {
              cat(class(object), "object\n")
              cat("Root:", object@root, "\n")
              cat("Number of bands: ", numNodes(object@toParentGraph), "\n")
              cat("Number of levels: ", length(object@level2nodes) - 1, "\n")
          })


setMethod("show", signature(object="LinearMResultBase"),
          function(object) {
              cat(description(object), "\n")
              nSig <- sum(pvalues(object) < object@pvalueCutoff[1])
              cat(length(pvalues(object)), testName(object), "ids tested ")
              cat("(", nSig, " have p < ", object@pvalueCutoff[1],
                  ")\n", sep="")
              cat("Selected gene set size:", geneMappedCount(object), "\n")
              cat("    Gene universe size:", universeMappedCount(object), "\n")
              cat("    Annotation package:", annotation(object), "\n")
          })


setMethod("show", signature(object="LinearMParams"),
          function(object) {
              cat("A", class(object), "instance\n")
              cat("  category:", object@categoryName, "\n")
              cat("annotation:", object@annotation, "\n")
          })
