## setMethod("pvalues", "ChrMapHyperGResult",
##           function(r) {
##               unlist(nodeData(r@chrGraph, attr="pvalue"))
##           })

## setMethod("universeMappedCount", "ChrMapHyperGResult",
##           function(r) {
##               u <- unlist(nodeData(r@chrGraph, attr="geneIds"))
##               length(unique(u))
##           })
setMethod("chrGraph", signature(r="ChrMapHyperGResult"),
          function(r) r@chrGraph)


setMethod("pvalues", signature(r="ChrMapHyperGResult"),
          function(r) {
              unlist(nodeData(r@chrGraph, attr="pvalue"))[r@pvalue.order]
              })


setMethod("oddsRatios", signature(r="ChrMapHyperGResult"),
          function(r) {
              unlist(nodeData(r@chrGraph, attr="oddsRatio"))[r@pvalue.order]
              })


setMethod("expectedCounts", signature(r="ChrMapHyperGResult"),
          function(r) {
              unlist(nodeData(r@chrGraph, attr="expCount"))[r@pvalue.order]
              })


entrezGeneUniverse <- function(r) {
    nodeData(r@chrGraph, n=nodes(r@chrGraph)[r@pvalue.order],
             attr="geneIds")
}


setMethod("geneIdUniverse", signature(r="ChrMapHyperGResult"),
          function(r) {
              entrezGeneUniverse(r)
          })


setMethod("condGeneIdUniverse", signature(r="ChrMapHyperGResult"),
          function(r) {
              if (isConditional(r))
                nodeData(r@chrGraph, n=nodes(r@chrGraph)[r@pvalue.order],
                         attr="condGeneIds")
              else
                geneIdUniverse(r)
          })


setMethod("geneCounts", signature(r="ChrMapHyperGResult"),
          function(r) {
              sapply(condGeneIdUniverse(r), function(x) {
                  sum(r@geneIds %in% x)
              })
          })


## setMethod("condGeneCounts", signature(r="ChrMapHyperGResult"),
##           function(r) {
##               sapply(condGeneIdUniverse(r), function(x) {
##                   sum(r@geneIds %in% x)
##               })
##           })


setMethod("universeCounts", signature(r="ChrMapHyperGResult"),
          function(r) {
              sapply(entrezGeneUniverse(r), length)
          })


setMethod("universeMappedCount", signature(r="ChrMapHyperGResult"),
          function(r) {
              length(unique(unlist(entrezGeneUniverse(r))))
          })


setMethod("isConditional", signature(r="ChrMapHyperGResult"),
          function(r) r@conditional)


setMethod("description", signature(object="ChrMapHyperGResult"),
          function(object) {
              cond <- "Conditional"
              if (!isConditional(object))
                cond <- ""
              desc <- paste("Gene to %s", cond, "test for %s-representation")
              desc <- sprintf(desc,
                              paste(testName(object), collapse=" "),
                              testDirection(object))
              desc
          })
