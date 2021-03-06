\name{cb_test}
\alias{cb_test}

\title{Chromosome Band Tree-Based Hypothesis Testing}
\description{
  \code{cb_test} is a flexible tool for discovering interesting
  chromosome bands relative to a selected gene list.  The function
  supports local and global tests which can be carried out in a top down
  or bottom up fashion on the tree of chromosome bands.
}
\usage{
cb_test(selids, chrtree, level, dir = c("up", "down"),
       type = c("local", "global"), next.pval = 0.05,
       cond.pval = 0.05, conditional = FALSE)
}

\arguments{
  
  \item{selids}{A vector of gene IDs.  The IDs should match those used
    to annotatate the \code{ChrBandTree} given by \code{chrtree}.  In most
    cases, these will be Entrez Gene IDs.}
  
  \item{chrtree}{A \code{ChrBandTree} object representing the chromosome
  bands and the mapping to gene identifiers.  The genes in the
  \code{ChrBandTree} are limited to the universe of gene IDs specified
  at object creation time.}

\item{level}{An integer giving the level of the chromosome band tree
  at which testing should begin.  The level is conceptualized as the set
  of nodes with a given path length to the root (organism) node of the
  chromosome band tree.  So level 1 is the chromosome and level 2 is the
  chromosome arms.  You can get a better sense by calling
  \code{exampleLevels(chrtree)}}

\item{dir}{A string giving the direction in which the chromosome band
  tree will be traversed when carrying out the tests.  A bottom up
  traversal, from leaves to root, is specified by \code{"up"}.  A top
  down, from root to leaves, traversal is specified by \code{"down"}.}

  \item{type}{A string giving the type of test to perform.  The current
    choices are \code{"local"} and \code{"global"}.  A local test
    carries out a chisq.test on each 2 x K contingency table induced by
    each set of siblings at a given level in the tree.  A global test
    uses the Hypergeometric distribution to compute a p-value for the 2
    x 2 tables induced by each band treated independently.}

  \item{next.pval}{The p-value cutoff used to determine whether the parents
    or children of a node should be tested.  After testing a given level
    of the tree, the decision of whether or not to continue testing the
    children (or parents) of the already tested nodes is made by
    comparing the p-value result for a given node with this cutoff;
    relatives of nodes with values strictly greater than the cutoff are
    skipped.}

  \item{cond.pval}{The p-value cutoff used to determine whether a node
    is significant during a conditional test.  See \code{conditional}.}

  \item{conditional}{A logical value.  Can only be used when
    \code{dir="up"} and \code{type="global"}.  In this case, a
    \code{TRUE} value causes a conditional Hypergeometric calculation to
    be performed.  The genes annotated at significant children of a
    given band are removed before testing.}
}

\value{
  A list with an element for each level of the tree that was tested.
  Note that the first element will correspond to the level given by
  \code{level} and that subsequent elements will be the next or previous
  depending on \code{dir}.

  Each level element is itself a list consisting of a result list for
  each node or set of nodes tested.  These inner-most lists will have,
  at least, the following components:

    \item{nodes}{A character vector of the nodes involved in the test.}
    \item{p.value}{The p-value for the test}
    \item{observed}{The contingency table}
    \item{method}{A brief description of the test method}
}

\author{Seth Falcon}


