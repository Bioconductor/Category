setMethod("show", signature(object="GeneCategoryHyperGeoTestResult"),
          function(object) {
              cat("Gene to", object@testName, 
                  "Category Association Test Result\n")
              nSig <- sum(object@pvalues < 0.05)
              cat(length(object@pvalues), object@testName, "ids tested ")
              cat("(", nSig, " have p < 0.05)\n", sep="")
              cat("Selected gene set size:", object@geneMappedCount, "\n")
              cat("    Gene universe size:", object@universeMappedCount, "\n")
              cat("    Annotation package:", object@annotation, "\n")
          })
