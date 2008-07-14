setMethod("geneIds", "ChrBandTree",
          function(object, n, ...) {
              stopifnot(length(n) == 1L)
              nodeData(object@toChildGraph, n=n, attr="geneIds")[[1L]]
          })

setMethod("lgeneIds", "ChrBandTree",
          function(r, n, ...) {
              ans <- nodeData(r@toChildGraph, n=n, attr="geneIds")
              ans[sapply(ans, function(x) !is.na(x[1]))]
          })

setMethod("childrenOf", signature=signature(g="ChrBandTree", n="character"),
          function(g, n, ...) {
              edges(g@toChildGraph)[n]
          })

setMethod("parentOf", signature=signature(g="ChrBandTree", n="character"),
          function(g, n, ...) {
              edges(g@toParentGraph)[n]
          })

setMethod("allGeneIds", "ChrBandTree",
          function(g, ...) {
              unique(unlist(nodeData(g@toChildGraph,
                                     attr="geneIds")))
          })

setMethod("treeLevels", "ChrBandTree",
          function(g, ...) {
              as.integer(names(g@level2nodes))
          })

setMethod("level2nodes", c("ChrBandTree", "numeric"),
          function(g, level, ...) {
              level2nodes(g, as.character(level), ...)
          })

setMethod("level2nodes", c("ChrBandTree", "character"),
          function(g, level, ...) {
              ans <- g@level2nodes[[level, exact=TRUE]]
              if (is.null(ans))
                stop("no such level: ", as.character(level))
              ans
          })

exampleLevels <- function(g)
{
    sapply(g@level2nodes, function(x) x[1])
}

make_level2node_map <- function(g)
{
    # Given a graph representing a tree with the root node given by
    # nodes(g)[1], return a list mapping levels of the tree to nodes.
    paths <- dijkstra.sp(g)
    split(names(paths[["distances"]]), paths[["distances"]])
}

ChrBandTreeFromGraph <- function(g) {
    root <- nodes(g)[1]
    gToParent <- reverseEdgeDirections(g)
    l2n <- make_level2node_map(g)
    new("ChrBandTree", toParentGraph=gToParent,
        toChildGraph=g, root=root, level2nodes=l2n)
}

NewChrBandTree <- function(chip, univ) {
    g <- makeChrBandGraph(chip, univ)
    ChrBandTreeFromGraph(g)
}
