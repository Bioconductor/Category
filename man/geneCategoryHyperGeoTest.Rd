\name{geneCategoryHyperGeoTest}
\docType{genericFunction}

\alias{geneCategoryHyperGeoTest}
\alias{geneCategoryHyperGeoTest,GeneCategoryHyperGeoTestParams-method}

\title{Test for overrepresentation of category terms in a set of genes}
\description{
  Given a subclass of \code{GeneCategoryHyperGeoTestParams}, compute
  Hypergeomtric p-values for overrepresentation of each term in the
  specified category among the specified gene set.
}
\usage{
geneCategoryHyperGeoTest(p)
}

\arguments{
  \item{p}{An instance of a subclass of
  \code{GeneCategoryHyperGeoTestParams}.  This parameter object
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
  '"YEAST"' then comparisons and statistics are computed using common
  names and are with respect to all genes annotated in the S. cerevisiae
  genome not with respect to any microarray chip.  This will *not* be
  the right thing to do if you are working with a yeast microarray.
}
\value{
  A \code{GeneCategoryHyperGeoTestResult} instance.
}

\section{Implementation Notes}{
  In most cases, the provided method with signature matching any
  subclass of \code{GeneCategoryHyperGeoTestParams} is all that will be
  needed.  This method follows a template pattern.  To add support for a
  new FOO category type, a developer would need to create a
  \code{GeneFooHyperGeoTestParams} subclass and then define two methods
  specialized to the new subclass that get called from inside
  \code{geneCategoryHyperGeoTest}: \code{universeBuilder} and
  \code{categoryToEntrezBuilder}.
}

\author{S. Falcon}
\seealso{
  \code{\link{GeneCategoryHyperGeoTestResult-class}}
  \code{\link{GeneCategoryHyperGeoTestParams-class}}
  \code{\link{GeneGoHyperGeoTestParams-class}}
  \code{\link{GeneKeggHyperGeoTestParams-class}}
}
\keyword{htest}
