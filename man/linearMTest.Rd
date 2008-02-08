\name{linearMTest}
\docType{genericFunction}

\alias{linearMTest}
\alias{linearMTest,LinearMParams-method}
% \alias{linearMTest,ChrMapLinearMParams-method}

\title{A linear model-based test to detect enrichment of unusual genes
  in categories}

\description{

  Given a subclass of \code{LinearMParams}, compute p-values for
  detecting systematic up or downregulation of the specified gene set in
  the specified category.

}

\usage{
  linearMTest(p)
}

\arguments{

  \item{p}{An instance of a subclass of \code{LinearMParams}.  This
    parameter object determines the category of interest (currently,
    only chromosome bands) as well as the gene set.
  }
}

\details{

  The per-gene statistics in the \code{geneStats} slot of \code{p} give
  a measure of up or down regulation of the individual genes in the
  universe.
%   The list of genes is reduced by removing identifiers that do not have
%   any annotations in the specified category.

%% FIXME: more details needed

%   It is important that the correct chip annotation data package be
%   identified as it determines the universe of gene identifiers and is
%   often used to determine the mapping between the category term and the
%   gene identifiers.
  
%   For S. cerevisiae if the \code{annotation} slot of \code{p} is set to
%   '"YEAST"' then comparisons and statistics are computed using common
%   names and are with respect to all genes annotated in the S. cerevisiae
%   genome not with respect to any microarray chip.  This will *not* be
%   the right thing to do if you are working with a yeast microarray.
}
\value{
  A \code{LinearMResult} instance.
}


\author{D. Sarkar}
\seealso{
  \code{\link{LinearMResult-class}}
  \code{\link{LinearMParams-class}}
}
\keyword{htest}
