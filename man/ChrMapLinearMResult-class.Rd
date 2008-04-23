\name{ChrMapLinearMResult-class}
\docType{class}
\alias{ChrMapLinearMResult-class}
\alias{chrGraph,ChrMapLinearMResult-method}
\alias{condGeneIdUniverse,ChrMapLinearMResult-method}
\alias{isConditional,ChrMapLinearMResult-method}

\title{Class "ChrMapLinearMResult" }
\description{

  This class represents the results of a linear model-based test for
  systematic changes in a per-gene statistic by chromosome band
  annotation.  The \code{\link{linearMTest}} function returns an
  instance of \code{ChrMapLinearMResult} when given a parameter object
  of class \code{ChrMapLinearMParams}.  Most slots can be queried using
  accessors.  }

\section{Objects from the Class}{

  Objects can be created by calls of the form
  \code{new("ChrMapLinearMResult", ...)}, but is more commonly created
  by callinf \code{\link{linearMTest}}

}

\section{Slots}{
  \describe{
    \item{\code{pvalue.order}:}{Object of class \code{"integer"} ~~ }
    \item{\code{conditional}:}{Object of class \code{"logical"} ~~ }
    \item{\code{chrGraph}:}{Object of class \code{"graph"} ~~ }
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
  Class \code{"\linkS4class{LinearMResult}"}, directly.

  Class \code{"\linkS4class{LinearMResultBase}"}, by class
  "LinearMResult", distance 2.
}
\section{Methods}{
  \describe{
    \item{chrGraph}{\code{signature(r = "ChrMapLinearMResult")}: ... }
    \item{condGeneIdUniverse}{\code{signature(r = "ChrMapLinearMResult")}: ... }
    \item{isConditional}{\code{signature(r = "ChrMapLinearMResult")}: ... }
	 }
}
\author{ Deepayan Sarkar }
\seealso{
  \code{\link{linearMTest}}, \code{\linkS4class{ChrMapLinearMParams}},
  \code{\linkS4class{LinearMResult}},
  \code{\linkS4class{LinearMResultBase}},
}
\examples{
showClass("ChrMapLinearMResult")
}
\keyword{classes}
