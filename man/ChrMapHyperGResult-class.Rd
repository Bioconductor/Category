\name{ChrMapHyperGResult-class}
\docType{class}
\alias{ChrMapHyperGResult-class}

\title{Class "ChrMapHyperGResult"}
\description{
  This class represents the results of a Hypergeometric test for
  over-representation of genes in a selected gene list in the
  chromosome band annotation.  The \code{hyperGTest} function returns
  an instance of \code{ChrMapHyperGResult} when given a parameter
  object of class \code{ChrMapHyperGParams}.  For details on accessing
  the results, see \link{HyperGResult-accessors}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ChrMapHyperGResult", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{pvalue.order}:}{Object of class \code{"integer"} that
    gives the order of the p-values.}
    \item{\code{conditional}:}{Object of class \code{"logical"}
    is a flag indicating whether or not this result is from a
    conditional analysis. }
    \item{\code{chrGraph}:}{Object of class \code{"graph"}.  The nodes
    are the chromosome bands with edges representing the tree
    structure of the bands.  Each node has a \code{"geneIds"}
    attribute that gives the gene IDs annotated at that band.}
    \item{\code{annotation}:}{A string giving the name of the chip
      annotation data package used.}
    \item{\code{geneIds}:}{Object of class \code{"ANY"}: the input
      vector of gene identifiers intersected with the universe of gene
      identifiers used in the computation.  The class of this slot is
      specified as \code{"ANY"} because gene IDs may be integer or
      character vectors depending on the annotation package.}
    \item{\code{testName}:}{A string identifying the testing method
      used.}
    \item{\code{pvalueCutoff}:}{Numeric value used a a p-value
        cutoff.  Used by the \code{show} method to count number of
        significant terms.
      }
    \item{\code{testDirection}:}{Object of class \code{"character"}
    indicating whether the test was for over-representation
    (\code{"over"}) or under-representation (\code{"under"}). }
  }
}
\section{Extends}{
Class \code{"\linkS4class{HyperGResultBase}"}, directly.
}
\section{Methods}{
See \link{HyperGResult-accessors}.
}

\author{Seth Falcon}
\examples{
showClass("ChrMapHyperGResult")
## For details on accessing the results:
##     help("HyperGResult-accessors")
}
\keyword{classes}
