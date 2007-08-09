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
          function(r, cond=TRUE) {
              if (cond && conditional(r))
                nodeData(r@chrGraph, n=nodes(r@chrGraph)[r@pvalue.order],
                         attr="condGeneIds")
              else
                entrezGeneUniverse(r)
          })


setMethod("condGeneIdUniverse", signature(r="ChrMapHyperGResult"),
          function(r) {
              geneIdUniverse(r, cond=TRUE)
          })


setMethod("isConditional", signature(r="ChrMapHyperGResult"),
          function(r) conditional(r))
