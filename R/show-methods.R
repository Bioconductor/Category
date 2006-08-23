setMethod("show", signature(object="GeneCategoryHyperGeoTestResult"),
          function(object) {
              cat("Gene to", object@testName, 
                  "Category Association Test Result\n")
              nSig <- sum(object@pvalues < object@pvalue.cutoff[1])
              cat(length(object@pvalues), object@testName, "ids tested ")
              cat("(", nSig, " have p < ", object@pvalue.cutoff[1],
                  ")\n", sep="")
              cat("Selected gene set size:", object@geneMappedCount, "\n")
              cat("    Gene universe size:", object@universeMappedCount, "\n")
              cat("    Annotation package:", object@annotation, "\n")
          })


setMethod("show", signature(object="GeneGoHyperGeoTestParams"),
          function(object) {
              cat("A", class(object), "instance\n")
              cat("  category:", object@categoryName, "\n")
              cat("annotation:", object@annotation, "\n")
          })
