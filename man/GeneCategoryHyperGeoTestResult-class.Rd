\name{GeneCategoryHyperGeoTestResult-class}
\docType{class}
\alias{GeneCategoryHyperGeoTestResult-class}
\alias{show,GeneCategoryHyperGeoTestResult-method}

\title{Class "GeneCategoryHyperGeoTestResult"}
\description{
  This class represents the results of a test for overrepresentation of
  categories among genes in a selected gene set based upon the
  Hypergeometric distribution.  The \code{geneCategoryHyperGeoTest}
  generic function returns an instance of the
  \code{GeneCategoryHyperGeoTestResult} class.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GeneCategoryHyperGeoTestResult", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{pvalues}:}{\code{"numeric"} vector: the ordered
      p-values for each category term tested.}
    \item{\code{geneCounts}:}{\code{"integer"} vector: for each
      category term tested, the number of genes from the gene set that
      are annotated at the term.}
    \item{\code{universeCounts}:}{\code{"integer"} vector: for
      each category term tested, the number of genes from the gene
      universe that are annotated at the term.}
    \item{\code{universeMappedCount}:}{\code{"integer"} vector of length
      one: The size of the gene universe set.}
    \item{\code{geneMappedCount}:}{\code{"integer"} vector of length
      one: The size of the selected gene set.}
    \item{\code{annotation}:}{A string giving the name of the chip
      annotation data package used.}
    \item{\code{geneIds}:}{Object of class \code{"ANY"}: the input
      vector of gene identifiers passed to the testing method.}
    \item{\code{testName}:}{A string identifying the testing method
      used.}
    \item{\code{pvalue.cutoff}:}{Numeric value used a a p-value
        cutoff.  Used by the \code{show} method to count number of
        significant terms.
    }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object =
    "GeneCategoryHyperGeoTestResult")}: Display a short summary of the result}
	 }
}

\author{S. Falcon}

\seealso{
  \code{\link{geneCategoryHyperGeoTest}}
  \code{\link{phyper}}
}

\keyword{classes}
