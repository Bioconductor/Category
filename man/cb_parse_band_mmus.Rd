\name{cb_parse_band_mmus}
\alias{cb_parse_band_mmus}

\title{Parse Mus Musculus Chromosome Band Annotations}
\description{
  This function parses chromosome band annotations as found in the
  <chip>MAP map of Bioconductor annotation data packages.  The return
  value is a vector of parent bands up to the relevant chromosome.
}
\usage{
cb_parse_band_mmus(x)
}

\arguments{
  \item{x}{A chromosome band annotation given as a string.}
}

\value{
  A character vector giving the path to the relevant chromosome.
}

\author{Seth Falcon}

\examples{
cb_parse_band_mmus("10 B3")
}

\keyword{ }

