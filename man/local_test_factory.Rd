\name{local_test_factory}
\alias{local_test_factory}
\alias{hg_test_factory}

\title{Local and Global Test Function Factories}
\description{
  These functions return functions appropriate for use as the
  \code{tfun} argument to \code{topdown_iter} or \code{bottomup_iter}.
  In particular, it is these functions that are associated with the
  "local" and "global" options for the \code{type} argument to
  \code{cb_test}.
}
\usage{
local_test_factory(selids, tableTest = chisq.test)
hg_test_factory(selids, PCUT = 0.05, COND = FALSE, OVER = TRUE)
}

\arguments{
  \item{selids}{A vector of gene IDs.  The IDs should match those used
    to annotatate the \code{ChrBandTree} given by \code{chrtree}.  In most
    cases, these will be Entrez Gene IDs.}
  \item{tableTest}{A contingency table testing function.  The behavior
    of this function must be reasonably close to that of \code{chisq.test}.}
  \item{PCUT}{A p-value cutoff that will be used to determine if a given
    test is significant or not when using \code{hg_test_factory}
    with \code{COND=TRUE}.}
  \item{COND}{A logical value indicating whether a conditional test
  should be performed.}
  \item{OVER}{If \code{TRUE}, test for over representation, if \code{FALSE},
    test for under representation.}
}
\details{
  The returned functions have signature \code{f(start, g, prev_ans)}
  where \code{start} is a vector of start nodes, \code{g} is a
  chromosome band tree graph, and \code{prev_ans} can contain the
  previous result returned by a call to this function.
}
\value{
  A function that can be used as the \code{tfun} argument to the tree
  visitor functions.
}

\author{Seth Falcon}
\seealso{\code{\link{cb_test}}}
\examples{
}

\keyword{internal}

