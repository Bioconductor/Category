\name{hyperGTest}
\docType{genericFunction}

\alias{hyperGTest}
\alias{hyperGTest,HyperGParams-method}
\alias{hyperGTest,ChrMapHyperGParams-method}
\alias{hyperGTest,KEGGHyperGParams-method}
\alias{hyperGTest,PFAMHyperGParams-method}

\title{Hypergeometric Test for association of categories and genes}

\description{
  Given a subclass of \code{HyperGParams}, compute Hypergeomtric
  p-values for over or under-representation of each term in the
  specified category among the specified gene set.
}
\usage{
hyperGTest(p)
}

\arguments{
  \item{p}{An instance of a subclass of
  \code{HyperGParams}.  This parameter object
  determines the category of interest (e.g., GO or KEGG) as well as the
  gene set.}
}
\details{
  The gene identifiers in the \code{geneIds} slot of \code{p} define the
  selected set of genes.  The universe of gene ids is determined by the
  chip annotation found in the \code{annotation} slot of \code{p}.  Both
  the selected genes and the universe are reduced by removing
  identifiers that do not have any annotations in the specified
  category.
  
  For each term in the specified category that has at least one
  annotation in the selected gene set, we determine how many of its
  annotations are in the universe set and how many are in the selected
  set.  With these counts we perform a Hypergeometric test using
  \code{phyper}.  This is equivalent to using Fisher's exact test.
  
  It is important that the correct chip annotation data package be
  identified as it determines the universe of gene identifiers and is
  often used to determine the mapping between the category term and the
  gene identifiers.
  
  For S. cerevisiae if the \code{annotation} slot of \code{p} is set to
  '"org.Sc.sgd"' then comparisons and statistics are computed using common
  names and are with respect to all genes annotated in the S. cerevisiae
  genome not with respect to any microarray chip.  This will *not* be
  the right thing to do if you are working with a yeast microarray.
}
\value{
  A \code{HyperGResult} instance.
}

\section{Implementation Notes}{
  In most cases, the provided method with signature matching any
  subclass of \code{HyperGParams} is all that will be needed.  This
  method follows a template pattern.  To add support for a new FOO
  category type, a developer would need to create a
  \code{FooHyperGParams} subclass and then define two methods
  specialized to the new subclass that get called from inside
  \code{hyperGTest}: \code{universeBuilder} and
  \code{categoryToEntrezBuilder}.
}

\author{S. Falcon}
\seealso{
  \code{\link{HyperGResult-class}}
  \code{\link{HyperGParams-class}}
  \code{\link{GOHyperGParams-class}}
  \code{\link{KEGGHyperGParams-class}}
}
\keyword{htest}
