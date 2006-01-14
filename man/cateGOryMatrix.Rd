\name{cateGOry}
\alias{cateGOry}
\title{Construct a category membership data structure from a list of gene identifiers and their
       annotated GO categories.}
\description{
  The function constructs a category membership data structure, such as used by
  \code{\link{applyByCategory}},
  from a list of gene identifiers and their annotated GO categories.
  For each of the GO categories stated in \code{categ}, it also
  considers the implied annotation with all ancestors.
}
\usage{
cateGOry(x, categ, sparse=TRUE)
}
\arguments{
  \item{x}{Character vector with (arbitrary) gene identifiers. They will
  be used for the column names of the resulting matrix.}
\item{categ}{A character vector of the same length as \code{x}
    with GO annotations for the genes in \code{x}. If a
    gene has multiple GO annotations, it is expected to occur multiple
    times in \code{x}, once for each different annotation.}
  \item{sparse}{Logical. See Section "Value".}
}
\details{
  Requires the \code{\link[GO:GO-package]{GO}} package.
}
\value{
  Depending on the value of \code{sparse}:
  
  If \code{TRUE}, an object of class
  \code{\link[graph:graphNEL-class]{graphNEL}}
  whose edges represent GO category membership of genes.
  This is really a bipartite graph, but the current implementation
  simply uses ordinary directed graph representation.

  If \code{FALSE}, an integer matrix whose columns correspond to
  \code{x} and whose rows correspond to GO categories. Row names are the
  GO category names, column names the gene identifiers. Matrix element
  \code{[i,j]} is 1 if gene \code{j} is in category \code{i},
  and 0 otherwise. For large vectors \code{x}, the sparse format will be
  more efficient.
}

\author{W. Huber}

\seealso{  \code{\link{applyByCategory}} }
\examples{
  g = cateGOry(c("CG2671", "CG2671", "CG2950"),
               c("GO:0000074", "GO:0001738", "GO:0003676"))
  edges(g)[c("GO:0000074",  "GO:0003676", "CG2950")]

  m = cateGOry(c("CG2671", "CG2671", "CG2950"),
               c("GO:0000074", "GO:0001738", "GO:0003676"), sparse=FALSE)

  for(ct in rownames(m))
    stopifnot(edges(g)[[ct]] == names(which(m[ct, ]!=0)))

}

\keyword{manip}

