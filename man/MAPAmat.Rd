\name{MAPAmat}
\alias{MAPAmat}
\title{ Create an adjacency matrix for MAP }
\description{
 For many organisms cytochrome bands, or MAP positions, can be determined
 on chromosomes. This function computes an adjacency matrix mapping 
 gene chip identifiers to MAP locations.
}
\usage{
MAPAmat(data, minCount = 5)
}
\arguments{
  \item{data}{ The character string naming the data package to use.}
  \item{minCount}{ The minimum number of genes that must be mapped
  to retain a position. }
}
\details{
 All MAP regions with fewer than \code{minCount} genes associated with it
 are dropped. Most genes are mapped to a single MAP region, in the rare
 cases where there are multiple mappings only the first one is used.
}
\value{
 A matrix of zeros and ones (possibly a sparse representation) where rows
 correspond to Entrez Gene IDs and columns to chromosome regions. A one 
  indicates that the corresponding gene is in the mapped region.

 Mapping is based on Entrez gene IDs.
}
\author{ R. Gentleman }
\seealso{ \code{\link{cateGOry}}, \code{\link{probes2MAP}} }
\examples{
have_hgu95av2 <- suppressWarnings(require("hgu95av2"))
if (have_hgu95av2)
   mam <- MAPAmat("hgu95av2", 10) 
}
\keyword{ manip }
