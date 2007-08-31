\name{tree_visitor}
\alias{tree_visitor}
\alias{bottomup_tree_visitor}
\alias{topdown_tree_visitor}

\title{Tree Visitor Function}
\description{
  This function visits each node in a tree-like object in an order
  determined by the \code{relationOf} function.  The function given by
  \code{tfun} is called for each set of nodes and the function
  \code{nfun} determines which nodes to test next optionally making use
  of the result of the previous test.
}
\usage{
tree_visitor(g, start, tfun, nfun, relationOf)
topdown_tree_visitor(g, start, tfun, nfun)
bottomup_tree_visitor(g, start, tfun, nfun)
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
    test.  This is used to determine the next set of nodes to test (see
    details).}

  \item{relationOf}{The method used to traverse the tree.  For example
    \code{childrenOf} or \code{parentOf}.}
}

\details{
  The \code{tree_visitor} function is intended to allow developers to
  quickly prototype different statistical testing paradigms on trees.
  It may be possible to extend this to work for DAGs.

  The visit begins by calling \code{tfun} with the nodes provided by
  \code{start}.  The result of each call to \code{tfun} is stored in an
  environment.  The concept is visitation by tree level and so each
  result is stored using a key representing the level (this isn't quite
  right since the nodes in \code{start} need not be first level, but
  they will be assigned key "1".  After storing the result, \code{nfun}
  is used to obtain a vector of accepted node labels.  The idea is that
  the user should have a way of determining which nodes in the next
  level of the tree are worth testing.  The next \code{start} set is
  determined by \code{start <- relationOf(g, accepted)} where accepted
  is \code{unique(nfun(ans, g))}.
}

\value{
  A list.  See the return value of \code{cb_test} to get an idea.  Each
  element of the list represents a call to \code{tfun} at a given level
  of the tree.
}
\author{Seth Falcon}

%% FIXME: if/when this gets exported, remove the internal keyword
\keyword{internal}

