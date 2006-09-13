\name{HyperGParams-class}
\docType{class}
\alias{HyperGParams-class}

\alias{categoryName}
\alias{categoryName,HyperGParams-method}

\title{Class "HyperGParams"}
\description{
  An abstract (VIRTUAL) parameter class for representing all parameters
  needed by a method specializing the \code{geneCategoryHyperGeoTest}
  generic.  You should only use subclasses of this class directly.
}
\section{Objects from the Class}{
Objects of this class cannot be instantiated directly.
}
\section{Slots}{
  \describe{
    \item{\code{geneIds}:}{Object of class \code{"ANY"}: A vector of
      gene identifiers.  Numeric and character vectors are probably the
      only things that make sense.  These are the gene ids for the
      selected gene set.}
    \item{\code{universeGeneIds}:}{Object of class \code{"ANY"}: A
      vector of gene ids in the same format as \code{geneIds} defining a
      subset of the gene ids on the chip that will be used as the
      universe for the hypergeometric calculation.  If this is
      \code{NULL} or has length zero, then all gene ids on the chip will
      be used.}
    \item{\code{annotation}:}{A string giving the name of the
      annotation data package for the chip used to generate the data.}
    \item{\code{cateogrySubsetIds}:}{Object of class \code{"ANY"}:
      If the test method supports it, can be used to specify a subset of
      category ids to include in the test instead of all possible
      category ids.}
    \item{\code{categoryName}:}{A string describing the category.
      Usually set automatically by subclasses.  For example "GO".}

    \item{\code{pvalueCutoff}:}{The p-value to use as a cutoff for
        significance for testing methods that require it.  This value
        will also be passed on to the result instance and used for
        display and counting of significant results.  The default is
        0.01.}

      \item{\code{testDirection}:}{A string indicating whether the
          test should be for overrepresentation (\code{"over"}) or
          underrepresentation (\code{"under"}).}
  }
}
\section{Methods}{
  \describe{
    \item{geneCategoryHyperGeoTest}{\code{signature(p =
    "HyperGParams")}: Perform hypergeometric tests to
    assess overrepresentation of category ids in the gene set.  See the
    documentation for the generic function for details.  This method
    must be called with a proper subclass of
    \code{HyperGParams}.}
}
}

\author{S. Falcon}

\seealso{
  \code{\link{HyperGResult-class}}
  \code{\link{GOHyperGParams-class}}
  \code{\link{KEGGHyperGParams-class}}
  \code{\link{geneKeggHyperGeoTest}}
  \code{\link{geneCategoryHyperGeoTest}}
}

\keyword{classes}
