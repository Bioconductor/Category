setMethod("show", signature(object="GeneCategoryHyperGeoTestResultBase"),
          function(object) {
              cat("Gene to", testName(object), 
                  "Category Association Test Result\n")
              cat("testDirection:", object@testDirection,
                  "representation\n")
              nSig <- sum(pvalues(object) < object@pvalueCutoff[1])
              cat(length(pvalues(object)), testName(object), "ids tested ")
              cat("(", nSig, " have p < ", object@pvalueCutoff[1],
                  ")\n", sep="")
              cat("Selected gene set size:", geneMappedCount(object), "\n")
              cat("    Gene universe size:", universeMappedCount(object), "\n")
              cat("    Annotation package:", annotation(object), "\n")
          })


setMethod("show", signature(object="GeneGoHyperGeoTestParams"),
          function(object) {
              cat("A", class(object), "instance\n")
              cat("  category:", object@categoryName, "\n")
              cat("annotation:", object@annotation, "\n")
          })
