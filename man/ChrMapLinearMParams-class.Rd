\name{ChrMapLinearMParams-class}
\docType{class}
\alias{ChrMapLinearMParams-class}
\alias{conditional<-,ChrMapLinearMParams,logical-method}
\alias{conditional,ChrMapLinearMParams-method}

\title{Class "ChrMapLinearMParams"}
\description{

  This class encapsulates parameters needed for testing systematic
  variations in some gene-level statistic by chromosome bands using
  \code{\link{linearMTest}}.

}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("ChrMapLinearMParams", ...)}.
}
\section{Slots}{
  \describe{

    \item{\code{chrGraph}:}{Object of class \code{"graph"}.  The nodes
    are the chromosome bands and the edges describe the tree structure
    of the bands.  Each node has a "geneIds" node attributes (see
    \code{nodeData}) which contains a vector of gene IDs annotated at
    the given band.}

    \item{\code{conditional}:}{Object of class \code{"logical"},
      indicating whether the test performed should be a conditional
      test.  }

    \item{\code{geneStats}:}{Named vector of class \code{"numeric"},
    giving the gene-level statistics to be used in the tests.  }
    
    \item{\code{universeGeneIds}:}{Object of class \code{"ANY"}: A
      vector of gene ids defining a subset of the gene ids on the chip
      that will be used as the universe for the hypergeometric
      calculation.  If this is \code{NULL} or has length zero, then all
      gene ids on the chip will be used.}

    \item{\code{annotation}:}{A string giving the name of the
      annotation data package for the chip used to generate the data.}

    \item{\code{datPkg}:}{Object of class \code{"DatPkg"} used to assist
      with dispatch based on type of annotation data available. }

    \item{\code{categorySubsetIds}:}{Object of class \code{"ANY"}:
      If the test method supports it, can be used to specify a subset of
      category ids to include in the test instead of all possible
      category ids.}

    \item{\code{categoryName}:}{A string describing the category.
      Usually set automatically by subclasses.  For example
      \code{"GO"}.}

    \item{\code{pvalueCutoff}:}{The p-value to use as a cutoff for
      significance for testing methods that require it.  This value will
      also be passed on to the result instance and used for display and
      counting of significant results.  The default is 0.01.}

    \item{\code{minSize}:}{An integer giving a minimum size for a gene
      set for it to be tested.  The default is 5. }

    \item{\code{testDirection}:}{A string indicating whether the test
      should test for systematic increase (\code{"up"}) or decrease
      (\code{"down"}) in the \code{geneStats} values within a gene set
      relative to the remaining genes.  }

  }
}
\section{Extends}{
  Class \code{"\linkS4class{LinearMParams}"}, directly.
}
\section{Methods}{
  \describe{

    \item{conditional}{\code{signature(r = "ChrMapLinearMParams")}:
      Accessor for the \code{conditional} slot. }

    \item{conditional<-}{\code{signature(r = "ChrMapLinearMParams",
	value = "logical")}: Replacement method for the
      \code{conditional} slot. }

  }
}
       
\author{ Deepayan Sarkar }

\seealso{
  \code{\link{linearMTest}}
}

\examples{
showClass("ChrMapLinearMParams")
}

\keyword{classes}
