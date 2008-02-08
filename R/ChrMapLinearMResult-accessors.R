setMethod("chrGraph", signature(r="ChrMapLinearMResult"),
          function(r) r@chrGraph)


setMethod("pvalues", signature(r="ChrMapLinearMResult"),
          function(r) {
              r@pvalues[r@pvalue.order]
          })


setMethod("geneIdUniverse", signature(r="ChrMapLinearMResult"),
          function(r, cond=TRUE) {
              r@catToGeneId[r@pvalue.order]
          })


setMethod("condGeneIdUniverse", signature(r="ChrMapLinearMResult"),
          function(r) {
              geneIdUniverse(r, cond=TRUE)
          })


setMethod("isConditional", signature(r="ChrMapLinearMResult"),
          function(r) conditional(r))
