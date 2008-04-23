\name{LinearMResult-class}
\docType{class}
\alias{LinearMResult-class}
\alias{effectSize,LinearMResult-method}
\alias{geneIdUniverse,LinearMResult-method}
\alias{pvalues,LinearMResult-method}

\title{Class "LinearMResult" }
\description{

  This class represents the results of a test for systematic change in
  some gene-level statistic by gene sets.  The \code{linearMTest}
  generic function returns an instance of the \code{LinearMResult}
  class. 

}
\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("LinearMResult",
    ...)}, but is more commonly created using a call to
    \code{\link{linearMTest}}.
  
}
\section{Slots}{
  \describe{
    \item{\code{pvalues}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{effectSize}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{catToGeneId}:}{Object of class \code{"list"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{geneIds}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{testName}:}{Object of class \code{"character"} ~~ }
    \item{\code{pvalueCutoff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{minSize}:}{Object of class \code{"integer"} ~~ }
    \item{\code{testDirection}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{LinearMResultBase}"}, directly.
}
\section{Methods}{
  \describe{
    \item{effectSize}{\code{signature(r = "LinearMResult")}: ... }
    \item{geneIdUniverse}{\code{signature(r = "LinearMResult")}: ... }
    \item{pvalues}{\code{signature(r = "LinearMResult")}: ... }
  }
}

\author{ Deepayan Sarkar }

\seealso{
  \code{\link{linearMTest}}
}
\examples{
showClass("LinearMResult")
}
\keyword{classes}
