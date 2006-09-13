\name{GOHyperGParams-class}
\docType{class}
\alias{GOHyperGParams-class}
\alias{show,GOHyperGParams-method}
\alias{categoryName,GOHyperGParams-method}

\title{Class "GOHyperGParams"}

\description{
  A parameter class for representing all parameters needed for running
  the \code{geneCategoryHyperGeoTest} method with one of the GO
  ontologies (BP, CC, MF) as the category.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GOHyperGParams", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{ontology}:}{A string specifying the GO ontology to use.
      Must be one of "BP", "CC", or "MF".}
    \item{\code{conditional}:}{A logical indicating whether the
      calculation should condition on the GO structure.}
    \item{\code{geneIds}:}{Object of class \code{"ANY"}: A vector of
      gene identifiers.  Numeric and character vectors are probably the
      only things that make sense.  These are the gene ids for the
      selected gene set.}
    \item{\code{universeGeneIds}:}{Object of class \code{"ANY"}: A
      vector of gene ids in the same format as \code{geneIds} defining a
      subset of the gene ids on the chip that will be used as the
      universe for the hypergeometric calculation.  If this is
      \code{NULL} or has length zero, then all gene ids on the chip will
      be used.}
    \item{\code{annotation}:}{A string giving the name of the
      annotation data package for the chip used to generate the data.}
    \item{\code{cateogrySubsetIds}:}{Object of class \code{"ANY"}:
      If the test method supports it, can be used to specify a subset of
      category ids to include in the test instead of all possible
      category ids.}
    \item{\code{categoryName}:}{A string describing the category.
      Usually set automatically by subclasses.  For example "GO".}           
  }
}
\section{Extends}{
Class \code{"HyperGParams"}, directly.
}
\section{Methods}{
  \describe{
    \item{geneCategoryHyperGeoTest}{\code{signature(p =
        "HyperGParams")}: Perform hypergeometric tests to
      assess overrepresentation of category ids in the gene set.  See
      the documentation for the generic function for details.  This
      method must be called with a proper subclass of
      \code{HyperGParams}.}
  }
}

\author{S. Falcon}

\seealso{
  \code{\link{HyperGResult-class}}
  \code{\link{GOHyperGParams-class}}
  \code{\link{geneKeggHyperGeoTest}}
  \code{\link{geneCategoryHyperGeoTest}}
}

\keyword{classes}
