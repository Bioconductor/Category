\name{LinearMResultBase-class}
\docType{class}
\alias{LinearMResultBase-class}
\alias{annotation,LinearMResultBase-method}
\alias{condGeneIdUniverse,LinearMResultBase-method}
\alias{conditional,LinearMResultBase-method}
\alias{description,LinearMResultBase-method}
\alias{geneIdsByCategory,LinearMResultBase-method}
\alias{geneIds,LinearMResultBase-method}
\alias{geneMappedCount,LinearMResultBase-method}
\alias{isConditional,LinearMResultBase-method}
\alias{pvalueCutoff,LinearMResultBase-method}
\alias{show,LinearMResultBase-method}
\alias{sigCategories,LinearMResultBase-method}
\alias{summary,LinearMResultBase-method}
\alias{testDirection,LinearMResultBase-method}
\alias{testName,LinearMResultBase-method}
\alias{universeCounts,LinearMResultBase-method}
\alias{universeMappedCount,LinearMResultBase-method}

\title{Class "LinearMResultBase" }
\description{	 
  This VIRTUAL class represents common elements of the return values
  of generic functions like \code{linearMTest}.
}

\section{Objects from the Class}{
  A virtual Class: No objects may be created from it.
}

\section{Slots}{
	 \describe{
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{geneIds}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{testName}:}{Object of class \code{"character"} ~~ }
    \item{\code{pvalueCutoff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{minSize}:}{Object of class \code{"integer"} ~~ }
    \item{\code{testDirection}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{annotation}{\code{signature(object = "LinearMResultBase")}: ... }
    \item{condGeneIdUniverse}{\code{signature(r = "LinearMResultBase")}: ... }
    \item{conditional}{\code{signature(r = "LinearMResultBase")}: ... }
    \item{description}{\code{signature(object = "LinearMResultBase")}: ... }
    \item{geneIdsByCategory}{\code{signature(object = "LinearMResultBase")}: ... }
    \item{geneIds}{\code{signature(object = "LinearMResultBase")}: ... }
    \item{geneMappedCount}{\code{signature(r = "LinearMResultBase")}: ... }
    \item{isConditional}{\code{signature(r = "LinearMResultBase")}: ... }
    \item{pvalueCutoff}{\code{signature(r = "LinearMResultBase")}: ... }
    \item{show}{\code{signature(object = "LinearMResultBase")}: ... }
    \item{sigCategories}{\code{signature(r = "LinearMResultBase")}: ... }
    \item{summary}{\code{signature(object = "LinearMResultBase")}: ... }
    \item{testDirection}{\code{signature(r = "LinearMResultBase")}: ... }
    \item{testName}{\code{signature(r = "LinearMResultBase")}: ... }
    \item{universeCounts}{\code{signature(r = "LinearMResultBase")}: ... }
    \item{universeMappedCount}{\code{signature(r = "LinearMResultBase")}: ... }
	 }
}
\author{ Deepayan Sarkar }

\seealso{
  \code{\linkS4class{LinearMResult}}, \code{\linkS4class{ChrMapLinearMParams}}
}
\examples{
showClass("LinearMResultBase")
}
\keyword{classes}
