\name{LinearMParams-class}
\docType{class}
\alias{LinearMParams-class}
\alias{annotation<-,LinearMParams,character-method}
\alias{annotation,LinearMParams-method}
\alias{categoryName,LinearMParams-method}
\alias{conditional,LinearMParams-method}
\alias{geneIds<-,LinearMParams-method}
\alias{geneIds,LinearMParams-method}
\alias{pvalueCutoff<-,LinearMParams-method}
\alias{pvalueCutoff,LinearMParams-method}
\alias{show,LinearMParams-method}
\alias{testDirection<-,LinearMParams-method}
\alias{testDirection,LinearMParams-method}
\alias{universeGeneIds,LinearMParams-method}

\title{Class "LinearMParams" }
\description{
  An abstract (VIRTUAL) parameter class for representing all parameters
  needed by a method specializing the \code{linearMTest} generic.  You
  should only use subclasses of this class directly.
}

\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{geneStats}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{universeGeneIds}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{datPkg}:}{Object of class \code{"DatPkg"} ~~ }
    \item{\code{cateogrySubsetIds}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{categoryName}:}{Object of class \code{"character"} ~~ }
    \item{\code{pvalueCutoff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{minSize}:}{Object of class \code{"integer"} ~~ }
    \item{\code{testDirection}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{annotation<-}{\code{signature(object = "LinearMParams", value = "character")}: ... }
    \item{annotation}{\code{signature(object = "LinearMParams")}: ... }
    \item{categoryName}{\code{signature(r = "LinearMParams")}: ... }
    \item{conditional}{\code{signature(r = "LinearMParams")}: ... }
    \item{geneIds<-}{\code{signature(object = "LinearMParams")}: ... }
    \item{geneIds}{\code{signature(object = "LinearMParams")}: ... }
    \item{pvalueCutoff<-}{\code{signature(r = "LinearMParams")}: ... }
    \item{pvalueCutoff}{\code{signature(r = "LinearMParams")}: ... }
    \item{show}{\code{signature(object = "LinearMParams")}: ... }
    \item{testDirection<-}{\code{signature(r = "LinearMParams")}: ... }
    \item{testDirection}{\code{signature(r = "LinearMParams")}: ... }
    \item{universeGeneIds}{\code{signature(r = "LinearMParams")}: ... }
  }
}
\author{ Deepayan Sarkar }

\seealso{

  \code{\linkS4class{ChrMapLinearMParams}} for descriptions of the
  slots.
  
}
\examples{
showClass("LinearMParams")
}
\keyword{classes}
