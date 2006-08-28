\name{GeneCategoryHyperGeoTestResultBase-class}
\docType{class}
\alias{GeneCategoryHyperGeoTestResultBase-class}

\alias{annotation,GeneCategoryHyperGeoTestResultBase-method}

\alias{geneIds}
\alias{geneIds,GeneCategoryHyperGeoTestResultBase-method}

\alias{geneMappedCount}
\alias{geneMappedCount,GeneCategoryHyperGeoTestResultBase-method}

\alias{pvalueCutoff}
\alias{pvalueCutoff,GeneCategoryHyperGeoTestResultBase-method}

\alias{show,GeneCategoryHyperGeoTestResultBase-method}

\alias{testName}
\alias{testName,GeneCategoryHyperGeoTestResultBase-method}

\title{Class "GeneCategoryHyperGeoTestResultBase"}
\description{	 
  This VIRTUAL class represents common elements of the return values
  of generic functions like \code{geneCategoryHyperGeoTest}.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
	 \describe{

         \item{\code{annotation}:}{Object of class \code{"character"}
             giving the name of the annotation data package used. }

         \item{\code{geneIds}:}{Object of class \code{"ANY"} (usually
             a character vector, but sometimes an integer vector).
             The input vector of gene identifiers intersected with the
             universe of gene identifiers used in the computation.  }

         \item{\code{testName}:}{Object of class \code{"character"}
             identifying the testing method used.
           }

         \item{\code{pvalue.cutoff}:}{Numeric value used by the
             testing method as a p-value cutoff.  Not all testing
             methods use this.  Also used by the \code{show} method to
             count number of significant terms.  }

  }
}
\section{Methods}{
  \describe{
    \item{annotation}{\code{signature(object =
          "GeneCategoryHyperGeoTestResultBase")}: return the name
        of the annotation data package used.  }

    \item{geneIds}{\code{signature(r =
          "GeneCategoryHyperGeoTestResultBase")}: return the input
        vector of gene identifiers intersected with the universe of
        gene identifiers used in the computation. }

    \item{geneMappedCount}{\code{signature(r =
          "GeneCategoryHyperGeoTestResultBase")}: return the size of
        the selected gene set used in the computation.  This is simply
        the \code{length(geneIds(obj))}.}

    \item{pvalueCutoff}{\code{signature(r =
          "GeneCategoryHyperGeoTestResultBase")}: accessor for the
        \code{pvalue.cutoff} slot. }

    \item{show}{\code{signature(object =
          "GeneCategoryHyperGeoTestResultBase")}: display an instance. }

    \item{testName}{\code{signature(r =
          "GeneCategoryHyperGeoTestResultBase")}: return a string
        describing the testing method used.}

	 }
}

\author{Seth Falcon}


\seealso{
  \code{\link{GeneCategoryHyperGeoTestResult-class}}
  \code{\link[GOstats]{GeneGoHyperGeoTestResult-class}}
}

\keyword{classes}
