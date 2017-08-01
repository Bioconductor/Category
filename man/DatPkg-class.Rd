\name{DatPkg-class}
\docType{class}
\alias{DatPkg-class}
\alias{AffyDatPkg-class}
\alias{YeastDatPkg-class}
\alias{Org.XX.egDatPkg-class}
\alias{GeneSetCollectionDatPkg}
\alias{OBOCollectionDatPkg-class}
\alias{OBOCollectionDatPkg}

\alias{DatPkgFactory}
\alias{DatPkgFactory,character-method}
\alias{DatPkgFactory,OBOCollection,GeneSetCollection-method}


\alias{ID2GO}
\alias{ID2KEGG}
\alias{ID2EntrezID}
\alias{GO2AllProbes}
\alias{KEGG2AllProbes}

\alias{ID2GO,DatPkg-method}
\alias{ID2GO,GeneSetCollectionDatPkg-method}
\alias{ID2KEGG,DatPkg-method}
\alias{ID2KEGG,GeneSetCollectionDatPkg-method}
\alias{ID2EntrezID,AffyDatPkg-method}
\alias{ID2EntrezID,YeastDatPkg-method}
\alias{ID2EntrezID,Org.XX.egDatPkg-method}
\alias{ID2EntrezID,ArabidopsisDatPkg-method}
\alias{ID2EntrezID,GeneSetCollectionDatPkg-method}
\alias{GO2AllProbes,DatPkg-method}
\alias{GO2AllProbes,YeastDatPkg-method}
\alias{GO2AllProbes,Org.XX.egDatPkg-method}
\alias{GO2AllProbes,GeneSetCollectionDatPkg-method}
\alias{KEGG2AllProbes,DatPkg-method}
\alias{KEGG2AllProbes,GeneSetCollectionDatPkg-method}
\alias{isDBDatPkg,DatPkg-method}
\alias{isDBDatPkg,GeneSetCollectionDatPkg-method}
\alias{organism,DatPkg-method}
\alias{organism,GeneSetCollectionDatPkg-method}


\title{Class "DatPkg"}
\description{
  \code{DatPkg} is a \code{VIRTUAL} class for representing annotation
  data packages.

  \code{AffyDatPkg} is a subclass of \code{DatPkg} used to represent
  standard annotation data packages that follow the format of Affymetrix
  expression array annotation.

  \code{YeastDatPkg} is a subclass of \code{DatPkg} used to represent
  the annotation data packages for yeast.  The yeast chip packages are
  based on sgd and are internally different from the \code{AffyDataPkg}
  conforming packages.

  \code{Org.XX.egDatPkg} is a subclass of \code{DatPkg} used to
  represent the \code{org.*.eg.db} organism-level Entez Gene based
  annotation data packages.
  
  \code{OBOCollectionDatPkg} is a subclass of \code{DatPkg} used to
  represent the \code{OBO} based annotation data packages.
  
  \code{GeneSetCollectionDatPkg} is a subclass of \code{DatPkg} used to
  represent annotations in the form of \code{GeneSetCollection} objects
  which are not based on any annotation packages but are instead derived
  from custom (user supplied) annotations.  
}
\section{Objects from the Class}{
  A virtual Class: No objects may be created from it.

  Given the name of an annotation data package, \code{DatPkgFactory} can
  be used to create an appropriate \code{DatPkg} subclass.
}
\section{Slots}{
  \describe{
    \item{name}{A string giving the name of the annotation data package.}
  }
}

\section{Methods}{
  See \code{showMethods(classes="DatPkg")}.
  %% FIXME: TODO: add doc for methods
  The set of methods, \code{ID2EntreizID} map between the standard
  IDs for an organism, or Chip and EntrezIDs, typically to give
  a way to get the GO terms.  Different organisms, such as S. cerevisae
  and A. thaliana have their own internal IDs, so we need specialized methods
  for them.
}

\author{Seth Falcon}
\examples{
DatPkgFactory("hgu95av2")
\dontrun{
DatPkgFactory("org.Sc.sgd")
DatPkgFactory("org.Hs.eg.db")
DatPkgFactory("ag")
}
}
\keyword{classes}
