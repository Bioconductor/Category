\name{HyperGParams-class}
\docType{class}

%% I don't know where these really belong :-(
\alias{conditional}
\alias{conditional<-}
\alias{geneIds<-}
\alias{ontology}
\alias{ontology<-}
\alias{pvalueCutoff<-}
\alias{testDirection<-}
\alias{universeGeneIds}
%% end :-(

\alias{HyperGParams-class}
\alias{initialize,HyperGParams-method}

\alias{geneIds,HyperGParams-method}
\alias{geneIds<-,HyperGParams,logical-method}
\alias{geneIds<-,HyperGParams,ANY-method}

\alias{annotation,HyperGParams-method}
\alias{annotation<-,HyperGParams,character-method}
\alias{conditional,HyperGParams-method}

\alias{ontology,HyperGParams-method}
\alias{pvalueCutoff,HyperGParams-method}
\alias{pvalueCutoff<-,HyperGParams-method}
\alias{testDirection,HyperGParams-method}
\alias{testDirection<-,HyperGParams-method}
\alias{universeGeneIds,HyperGParams-method}

\alias{categoryName}
\alias{categoryName,HyperGParams-method}

\alias{organism,HyperGParams-method}

\alias{makeValidParams,HyperGParams-method}
\alias{show,HyperGParams-method}

\title{Class "HyperGParams"}
\description{
  An abstract (VIRTUAL) parameter class for representing all parameters
  needed by a method specializing the \code{hyperGTest}
  generic.  You should only use subclasses of this class directly.
}
\section{Objects from the Class}{
Objects of this class cannot be instantiated directly.
}
\section{Slots}{
  \describe{
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
    \item{\code{annotation}:}{Object of class
      \code{"ANY"}. Functionally, this is either a string giving the name of the
      annotation data package for the chip used to generate the data, or
      the name of an annotation object downloaded using AnnotationHub.}
    \item{\code{categorySubsetIds}:}{Object of class \code{"ANY"}:
      If the test method supports it, can be used to specify a subset of
      category ids to include in the test instead of all possible
      category ids.}
    \item{\code{categoryName}:}{A string describing the category.
      Usually set automatically by subclasses.  For example "GO".}
    \item{\code{pvalueCutoff}:}{The p-value to use as a cutoff for
      significance for testing methods that require it.  This value
      will also be passed on to the result instance and used for
      display and counting of significant results.  The default is
      0.01.}
    \item{\code{testDirection}:}{A string indicating whether the test
      should be for overrepresentation (\code{"over"}) or
      underrepresentation (\code{"under"}).}
    \item{\code{datPkg}:}{Holds a DatPkg object which is of a
     particular type that in turn varies with the kind of annotation
     package this is.}

  }
}
\section{Methods}{
  \describe{
  \item{hyperGTest}{\code{signature(p =
        "HyperGParams")}: Perform hypergeometric tests to
      assess overrepresentation of category ids in the gene set.  See the
      documentation for the generic function for details.  This method
      must be called with a proper subclass of
      \code{HyperGParams}.}
    
  \item{\code{geneIds(object)}, \code{geneIds(object) <- value}}{Accessors for
      the gene identifiers that will be used as the selected gene
      list.}
    
  \item{code{annotation(object)}}{Accessor for annotation.  If you want
    to change the annotation for an existing instance, use the
    replacement form.}
    
  \item{\code{ontology(object)}}{Accessor for GO ontology.}

  \item{\code{organism(object)}}{Accessor for the organism character
    string used as an identifier in \code{DatPkg}.}
    
  \item{\code{pvalueCutoff(r)}, \code{pvalueCutoff(r) <-
        value}}{Accessor for the p-value cutoff.  When setting,
      \code{value} should be a numeric value between zero and one.}

  \item{\code{testDirection}}{Accessor for the test direction.  When setting,
        \code{value} must be either "over" or "under".}

    \item{\code{universeGeneIds(r)}}{accessor for vector of gene
        identifiers.}

  }
}

\author{S. Falcon}

\seealso{
  \code{\link{HyperGResult-class}}
  \code{\link{GOHyperGParams-class}}
  \code{\link{KEGGHyperGParams-class}}
  \code{\link{hyperGTest}}
}

\keyword{classes}
