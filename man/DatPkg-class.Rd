\name{DatPkg-class}
\docType{class}
\alias{DatPkg-class}
\alias{AffyDatPkg-class}
\alias{YeastDatPkg-class}
\alias{OrganismMappingDatPkg-class}
\alias{Org.XX.egDatPkg-class}

\alias{DatPkgFactory}

\alias{ID2GO}
\alias{ID2EntrezID}
\alias{GO2AllProbes}

\alias{ID2GO,DatPkg-method}
\alias{ID2GO,OrganismMappingDatPkg-method}
\alias{ID2EntrezID,AffyDatPkg-method}
\alias{ID2EntrezID,YeastDatPkg-method}
\alias{ID2EntrezID,OrganismMappingDatPkg-method}
\alias{ID2EntrezID,Org.XX.egDatPkg-method}
\alias{GO2AllProbes,DatPkg-method}
\alias{GO2AllProbes,OrganismMappingDatPkg-method}
\alias{GO2AllProbes,Org.XX.egDatPkg-method}

\title{Class "DatPkg"}
\description{
  \code{DatPkg} is a \code{VIRTUAL} class for representing annotation
  data packages.

  \code{AffyDatPkg} is a subclass of \code{DatPkg} used to represent
  standard annotation data packages that follow the format of Affymetrix
  expression array annotation.

  \code{YeastDatPkg} is a subclass of \code{DatPkg} used to represent
  the \code{YEAST} annotation data package.  The \code{YEAST} package
  needs a class because its format differs slightly from the format of
  \code{AffyDataPkg} conforming packages.

  \code{OrganismMappingDatPkg} is a subclass of \code{DatPkg} used to
  represent the \code{*LLMappings} annotation data packages.  The
  format of the \code{*LLMappings} packages is unique and this class
  allows for interface methods to make it interchangeable with the other
  formats.
  
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
}

\author{Seth Falcon}
\examples{
DatPkgFactory("hgu95av2.db")
DatPkgFactory("YEAST")
DatPkgFactory("fooLLMappings")

}
\keyword{classes}
