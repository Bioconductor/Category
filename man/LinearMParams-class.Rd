\name{LinearMParams-class}
\docType{class}
\alias{LinearMParams-class}
\alias{annotation<-,LinearMParams,character-method}
\alias{annotation,LinearMParams-method}
\alias{categoryName,LinearMParams-method}
\alias{conditional,LinearMParams-method}
\alias{conditional<-,LinearMParams,logical-method}
\alias{geneIds<-,LinearMParams,ANY-method}
\alias{geneIds,LinearMParams-method}
\alias{pvalueCutoff<-,LinearMParams-method}
\alias{pvalueCutoff,LinearMParams-method}
\alias{show,LinearMParams-method}
\alias{testDirection<-,LinearMParams-method}
\alias{testDirection,LinearMParams-method}
\alias{universeGeneIds,LinearMParams-method}

\title{Class "LinearMParams" }
\description{
  A parameter class for representing all parameters
  needed by a method specializing the \code{\link{linearMTest}} generic.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("LinearMParams", ...)}.}

\section{Slots}{
  \describe{
    \item{\code{geneStats}:}{Named vector of class \code{"numeric"},
      giving the gene-level statistics to be used in the tests. The
      names should correspond to the gene identifiers in \code{gsc}. }

    \item{\code{universeGeneIds}:}{Object of class \code{"ANY"}: A
      vector of gene ids defining a subset of the gene ids on the chip
      that will be used as the universe for the hypergeometric
      calculation.  If this is \code{NULL} or has length zero, then all
      gene ids on the chip will be used. Currently this parameter is
      ignored by the base \code{linearMTest} method.}

    \item{\code{annotation}:}{A string giving the name of the
      annotation data package for the chip used to generate the data.}
    
    \item{\code{datPkg}:}{Object of class \code{"DatPkg"} used to assist
      with dispatch based on type of annotation data
      available. Currently this parameter is ignored by the
      base \code{linearMTest} method.}
    
    \item{\code{categorySubsetIds}:}{Object of class \code{"ANY"}:
      If the test method supports it, can be used to specify a subset of
      category ids to include in the test instead of all possible
      category ids. Currently this parameter is ignored by the
      base \code{linearMTest} method.}

    \item{\code{categoryName}:}{A string describing the category.
      Usually set automatically by subclasses.  For example
      \code{"ChrMap"}. }

    \item{\code{pvalueCutoff}:}{The p-value to use as a cutoff for
      significance for testing methods that require it.  This value will
      also be passed on to the result instance and used for display and
      counting of significant results.  The default is 0.01. }

    \item{\code{minSize}:}{An integer giving a minimum size for a gene
      set for it to be tested.  The default is 5.  }

    \item{\code{testDirection}:}{A string indicating whether the test
      should test for systematic increase (\code{"up"}) or decrease
      (\code{"down"}) in the \code{geneStats} values within a gene set
      relative to the remaining genes.  }

    \item{\code{graph}:}{The \code{\link[graph:graph-class]{graph}}
      object indicating the hierarchical relationship among terms of the
      ontology or other grouping.
    }

    \item{\code{conditional}:}{A \code{logical} indicating whether
      conditional tests should be performed. This tests whether a term
      is still significant even when including its sub-terms in the model.
    }

    \item{\code{gsc}:}{The
      \code{\link[GSEABase:GeneSetCollection-class]{GeneSetCollection}}
      object grouping the gene ids into sets.
    }
  }
}
\section{Methods}{
  These are accessor methods for the various parameter slots:
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
    \item{conditional<-}{\code{signature(r = "LinearMParams")}: ... }
    \item{conditional}{\code{signature(r = "LinearMParams")}: ... }
    \item{universeGeneIds}{\code{signature(r = "LinearMParams")}: ... }
  }
}
\author{ Deepayan Sarkar, Michael Lawrence }

\seealso{
  See \code{\link{linearMTest}} for
  examples. \code{\linkS4class{ChrMapLinearMParams}} is a specialization
  of this class for chromosome maps.
}

\keyword{classes}
