\name{cb_parse_band_hsa}
\alias{cb_parse_band_hsa}

\title{Parse Homo Sapiens Chromosome Band Annotations}
\description{
  This function parses chromosome band annotations as found in the
  <chip>MAP map of Bioconductor annotation data packages.  The return
  value is a vector of parent bands up to the relevant chromosome.
}
\usage{
cb_parse_band_hsa(x)
}

\arguments{
  \item{x}{A chromosome band annotation given as a string.}
}

\value{
  A character vector giving the path to the relevant chromosome.
}

\author{Seth Falcon}

\examples{
cb_parse_band_hsa("12q32.12")
}

\keyword{ }

