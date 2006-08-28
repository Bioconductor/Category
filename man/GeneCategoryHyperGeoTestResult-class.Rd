\name{GeneCategoryHyperGeoTestResult-class}
\docType{class}
\alias{GeneCategoryHyperGeoTestResult-class}

\alias{geneCounts}
\alias{geneCounts,GeneCategoryHyperGeoTestResult-method}

\alias{pvalues}
\alias{pvalues,GeneCategoryHyperGeoTestResult-method}

\alias{universeCounts}
\alias{universeCounts,GeneCategoryHyperGeoTestResult-method}

\alias{universeMappedCount}
\alias{universeMappedCount,GeneCategoryHyperGeoTestResult-method}

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
    \item{\code{pvalue.cutoff}:}{Numeric value used a a p-value
        cutoff.  Used by the \code{show} method to count number of
        significant terms.
      }
    }
  }

\section{Extends}{
Class \code{"GeneCategoryHyperGeoTestResultBase"}, directly.
}
\section{Methods}{
  \describe{

  \item{geneCounts}{\code{signature(r =
        "GeneCategoryHyperGeoTestResult")}: return an \code{"integer"}
      vector: for each category term tested, the number of genes from
      the gene set that are annotated at the term.}

  \item{pvalues}{\code{signature(r =
        "GeneCategoryHyperGeoTestResult")}: return a \code{"numeric"}
      vector: the ordered p-values for each category term tested.}

  \item{universeCounts}{\code{signature(r =
        "GeneCategoryHyperGeoTestResult")}: return an \code{"integer"}
      vector: for each category term tested, the number of genes from
      the gene universe that are annotated at the term.}

  \item{universeMappedCount}{\code{signature(r =
        "GeneCategoryHyperGeoTestResult")}: return an \code{"integer"}
      vector of length one giving the size of the gene universe set. }


  }
}

\author{Seth Falcon}

\seealso{
  \code{\link{GeneCategoryHyperGeoTestResultBase-class}}
  \code{\link[GOstats]{GeneGoHyperGeoTestResult-class}}
}

\keyword{classes}
