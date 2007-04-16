\name{HyperGResultBase-class}
\docType{class}
\alias{HyperGResultBase-class}


\alias{show,HyperGResultBase-method}


\title{Class "HyperGResultBase"}
\description{	 
  This VIRTUAL class represents common elements of the return values
  of generic functions like \code{hyperGTest}.  All subclasses are
  intended to implement the accessor functions documented at
  \link{HyperGResult-accessors}.
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
  See \link{HyperGResult-accessors}.
}

\author{Seth Falcon}


\seealso{
  \code{\link{HyperGResult-class}}
  \code{\link[GOstats]{GeneGoHyperGeoTestResult-class}}
  \link{HyperGResult-accessors}
}

\keyword{classes}
