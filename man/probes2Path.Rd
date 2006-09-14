\name{probes2Path}
\alias{probes2Path}
\title{A function to map probe identifiers to pathways. }
\description{
 Given a set of probe identifiers from a microarray this function looks
 up all KEGG pathways that the probe is documented to be involved in.
}
\usage{
probes2Path(pids, data = "hgu133plus2")
}

\arguments{
  \item{pids}{A vector of probe identifiers. }
  \item{data}{The character name of the chip. }
}
\details{
  This is a simple look up in the appropriate chip \code{PATH} data
  environment. 
}
\value{
 A list of pathway vectors. One element for each value of \code{pid}
 that is mapped to at least one pathway.
}
\author{R. Gentleman }

\seealso{ \code{\link{findAMstats}} }
\examples{
 library("hgu95av2")
 x = c("1001_at", "1000_at")
 probes2Path(x, "hgu95av2")

}
\keyword{manip}
