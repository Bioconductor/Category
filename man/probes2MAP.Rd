\name{probes2MAP}
\alias{probes2MAP}
\title{ Map probe IDs to MAP regions. }
\description{
 This function maps probe identifiers to MAP positions using
 the appropriate Bioconductor meta-data package.
}
\usage{
probes2MAP(pids, data = "hgu133plus2")
}
\arguments{
  \item{pids}{ A vector of probe IDs for the chip in use.}
  \item{data}{ The name of the chip, as a character string. }
}
\details{
 Probes are mapped to regions, no checking for duplicate Entrez gene IDs 
 is done.
}
\value{
 A vector, the same length as \code{pids}, with the MAP locations.
}
\author{ R. Gentleman }
\seealso{ \code{\link{probes2Path}} }
\examples{
  set.seed(123)
  library("hgu95av2")
  v1 = sample(names(as.list(hgu95av2MAP)), 100)
  pp = probes2MAP(v1, "hgu95av2")
}
\keyword{ manip }
