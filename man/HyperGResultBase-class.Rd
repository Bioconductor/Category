\name{HyperGResultBase-class}
\docType{class}
\alias{HyperGResultBase-class}

\alias{annotation,HyperGResultBase-method}

\alias{geneIds}
\alias{geneIds,HyperGResultBase-method}

\alias{geneMappedCount}
\alias{geneMappedCount,HyperGResultBase-method}

\alias{pvalueCutoff}
\alias{pvalueCutoff,HyperGResultBase-method}

\alias{testDirection}
\alias{testDirection,HyperGResultBase-method}

\alias{show,HyperGResultBase-method}

\alias{description,HyperGResultBase-method}

\alias{testName}
\alias{testName,HyperGResultBase-method}

\title{Class "HyperGResultBase"}
\description{	 
  This VIRTUAL class represents common elements of the return values
  of generic functions like \code{hyperGTest}.
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

         \item{\code{pvalueCutoff}:}{Numeric value used by the
             testing method as a p-value cutoff.  Not all testing
             methods use this.  Also used by the \code{show} method to
             count number of significant terms.  }

         \item{\code{testDirection}:}{A string indicating whether the test
             performed was for overrepresentation (\code{"over"}) or
             underrepresentation(\code{"under"}).}

  }
}
\section{Methods}{
  \describe{
    \item{annotation}{\code{signature(object =
          "HyperGResultBase")}: return the name
        of the annotation data package used.  }

    \item{geneIds}{\code{signature(r =
          "HyperGResultBase")}: return the input
        vector of gene identifiers intersected with the universe of
        gene identifiers used in the computation. }

    \item{geneMappedCount}{\code{signature(r =
          "HyperGResultBase")}: return the size of
        the selected gene set used in the computation.  This is simply
        the \code{length(geneIds(obj))}.}

    \item{pvalueCutoff}{\code{signature(r =
          "HyperGResultBase")}: accessor for the
        \code{pvalueCutoff} slot. }

    \item{testDirection}{\code{signature(r =
          "HyperGResultBase")}: accessor for the
        \code{testDirection} slot.  Contains a string indicating
        whether the test was for \code{"over"} or \code{"under"}
        representation of the categories.}

    \item{description}{\code{signature(object =
          "HyperGResultBase")}: return a character
        string description of the test result. }

    \item{show}{\code{signature(object =
          "HyperGResultBase")}: display an instance. }

    \item{testName}{\code{signature(r =
          "HyperGResultBase")}: return a string
        describing the testing method used.}

	 }
}

\author{Seth Falcon}


\seealso{
  \code{\link{HyperGResult-class}}
  \code{\link[GOstats]{GeneGoHyperGeoTestResult-class}}
}

\keyword{classes}
