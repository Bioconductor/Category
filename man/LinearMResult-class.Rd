\name{LinearMResult-class}
\docType{class}
\alias{LinearMResult-class}
\alias{effectSize,LinearMResult-method}
\alias{summary,LinearMResult-method}
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
    \item{\code{pvalues}:}{Object of class \code{"numeric"}, with the
      p-values for each term.
    }
    \item{\code{pvalue.order}:}{Object of class \code{"integer"}, the
      order vector (increasing) for the p-values.
    }
    \item{\code{effectSize}:}{Object of class \code{"numeric"}, with
      the effect size for each term.
    }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{geneIds}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{testName}:}{Object of class \code{"character"} ~~ }
    \item{\code{pvalueCutoff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{minSize}:}{Object of class \code{"integer"} ~~ }
    \item{\code{testDirection}:}{Object of class \code{"character"} ~~ }
    \item{\code{conditional}:}{Object of class \code{"logical"} ~~ }
    \item{\code{graph}:}{Object of class \code{"graph"} ~~ }
    \item{\code{gsc}:}{Object of class \code{"GeneSetCollection"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{LinearMResultBase}"}, directly.
}
\section{Methods}{
  \describe{
    \item{effectSize}{\code{signature(r = "LinearMResult")}: ... }
    \item{pvalues}{\code{signature(r = "LinearMResult")}: ... }
    \item{summary}{\code{signature(r = "LinearMResult")}: returns
      a \code{data.frame} with a row for each gene set tested the
      following columns: \code{ID}, \code{Pvalue}, \code{Effect} size, 
      \code{Size} (number of members), \code{Conditional} (whether the
      test used the conditional test), and \code{TestDirection} (for up
      or down).
    }
  }
}

\author{ Deepayan Sarkar, Michael Lawrence }

\seealso{
  \code{\link{linearMTest}}
}
\examples{
showClass("LinearMResult")
}
\keyword{classes}
