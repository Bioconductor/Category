\name{HyperGResult-class}
\docType{class}
\alias{HyperGResult-class}

\title{Class "HyperGResult"}
\description{
  This class represents the results of a test for over-representation of
  categories among genes in a selected gene set based upon the
  Hypergeometric distribution.  The \code{hyperGTest}
  generic function returns an instance of the
  \code{HyperGResult} class. For details on accessing
  the results, see \link{HyperGResult-accessors}.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HyperGResult", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{pvalues}:}{\code{"numeric"} vector: the ordered
      p-values for each category term tested.}
    \item{\code{catToGeneId}:}{Object of class \code{"list"}.  The
        names of the list are category IDs.  Each element is a vector
        of gene IDs annotated at the given category ID and in the
        specified gene universe.}
    \item{\code{annotation}:}{A string giving the name of the chip
      annotation data package used.}
    \item{\code{geneIds}:}{Object of class \code{"ANY"}: the input
      vector of gene identifiers intersected with the universe of gene
      identifiers used in the computation.  The class of this slot is
      specified as \code{"ANY"} because gene IDs may be integer or
      character vectors depending on the annotation package.}
    \item{\code{testName}:}{A string identifying the testing method
      used.}
    \item{\code{pvalueCutoff}:}{Numeric value used a a p-value
        cutoff.  Used by the \code{show} method to count number of
        significant terms.
      }

      \item{\code{testDirection}:}{A string indicating whether the
          test should be for overrepresentation (\code{"over"}) or
          underrepresentation (\code{"under"}).}

  \item{\code{oddsRatios}}{a \code{"numeric"} vector giving
      the odds ratio for each category term tested.}

  \item{expectedCounts}{a \code{"numeric"} vector
      giving the expected number of genes in the selected gene list to
      be found at each tested category term.}


    }
  }

\section{Extends}{
Class \code{"HyperGResultBase"}, directly.
}
\section{Methods}{
See \link{HyperGResult-accessors}.
}

\author{Seth Falcon}

\seealso{
  \code{\link{HyperGResultBase-class}}
  \code{\link[GOstats]{GeneGoHyperGeoTestResult-class}}
  \link{HyperGResult-accessors}
}

\keyword{classes}
