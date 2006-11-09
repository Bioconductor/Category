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

\alias{isConditional,HyperGResultBase-method}

\alias{summary,HyperGResultBase-method}

\alias{htmlReport}
\alias{htmlReport,HyperGResultBase-method}

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

    \item{description}{\code{signature(object = "HyperGResultBase")}: return a character
        string description of the test result. }

    \item{show}{\code{signature(object =
          "HyperGResultBase")}: display an instance. }

    \item{testName}{\code{signature(r =
          "HyperGResultBase")}: return a string
        describing the testing method used.}

      \item{\code{isConditional(r)}}{Return \code{TRUE} if the
          result was obtained using a conditional algorithm.}

      \item{summary}{\code{signature(r = "HyperGResultBase")}: Returns
          a \code{data.frame} summarizing the test result.  Optional
          arguments \code{pvalue} and \code{categorySize} allow
          specification of maximum p-value and minimum cateogyrSize,
          respectively.}

      \item{htmlReport}{\code{signature(r = "HyperGResultBase")}:
          Write an HTML version of the table produced by the
          \code{summary} method.  The path of a file to write the
          report to can be specified using the \code{file} argument.
          The default is \code{file=""} which will cause the report to
          be printed to the screen.  If you wish to create a single
          report comprising multiple results you can set
          \code{append=TRUE}.  The default is \code{FALSE} (overwrite
          preexisting report file).  You can specify a string to use
          as an identifier for each table by providing a value for the
          \code{label} argument.  Additional named arguments will be
          passed to the \code{summary} method.}

    }
}

\author{Seth Falcon}


\seealso{
  \code{\link{HyperGResult-class}}
  \code{\link[GOstats]{GeneGoHyperGeoTestResult-class}}
}

\keyword{classes}
