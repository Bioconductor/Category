\name{effectSize}
\alias{effectSize}
%% move to a LinearMResult-accessors.Rd?
\docType{genericFunction}
\title{ Extract estimated effect sizes }
\description{
  This function extracts estimated effect sizes from the results of a
  linear model-based gene-set / category enrichment test. 
}
\usage{
effectSize(r)
}
\arguments{
  \item{r}{ The results of the test }
}
\value{
  A numeric vector.
}
\author{ Deepayan Sarkar }
\seealso{ \code{linkS4class{LinearMResult}} }
