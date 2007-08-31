\name{tree_iter}
\alias{tree_iter}
\alias{bottomup_iter}
\alias{topdown_iter}

\title{Tree Visitor Function}
\description{
  This function visits each node in a tree-like object in an order
  determined by the \code{relationOf} function.
}
\usage{
tree_iter(g, start, tfun, nfun, relationOf)
}

\arguments{

  \item{g}{A tree-like object that supports the method given by
  \code{relationOf}.}

  \item{start}{The set of nodes to start the computation (can be a list
    of siblings), but the nodes should all belong to the same level of
    the tree (same path length to root node).}

  \item{tfun}{The test function applied to each list of siblings at each
    level starting with \code{start}.  The signature of \code{tfun}
    should be \code{(start, g, prev_ans)}.}

  \item{nfun}{A function with signature \code{(ans, g)} that processes
    the result of \code{tfun} and returns a character vector of node
    names corresponding to nodes that were involved in an "interesting"
    test.  This is used to determine which }

  \item{relationOf}{The method used to traverse the tree.  For example
    \code{childrenOf} or \code{parentOf}.}
}

\value{
  See the return value of \code{cb_test} to get an idea.
}
\author{Seth Falcon}

\keyword{ }

