\name{getPathNames}
\alias{getPathNames}
\title{A function to print pathway names given their numeric id.}
\description{
 Given a KEGG pathway ID this function returns the character name of the
 pathway. 
}
\usage{
getPathNames(iPW)
}
\arguments{
  \item{iPW}{A vector of KEGG pathway IDs. }
}
\details{
  This function simply does a look up in \code{KEGGPATHID2NAME} and
  returns a list of the pathway names.

  Possible extensions would be to extend it to work with the cMAP
  library as well.
}
\value{
  A list of pathway names.
}

\author{R. Gentleman }

\seealso{\code{\link[KEGG]{KEGGPATHID2NAME}} }
\examples{

 nms = "00031"
 getPathNames(nms) 
}
\keyword{manip}

