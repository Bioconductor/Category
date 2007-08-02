\name{exampleLevels}
\alias{exampleLevels}

\title{Display a sample node from each level of a ChrBandTree object}
\description{
  The "levels" of a chromosome band tree represented by a \code{ChrBandTree} object
  are the sets of nodes with a given path length to the root node.  This
  function displays the available levels along with an example node from
  each level.
}
\usage{
exampleLevels(g)
}

\arguments{
  \item{g}{A \code{ChrBandTree} object}
}
\value{
  A list with an element for each level.  The names of the list are the
  levels.  Each element is an example of a node from the given level.
}

\author{S. Falcon}
% R documentation directory.
\keyword{ }

