\name{cateGOryMatrix}
\alias{cateGOryMatrix}
\title{Construct a category incidence matrix from a list of gene identifiers and their
       annotated GO categories.}
\description{
  The function constructs a category incidence matrix, such as used by
  \code{\link{applyByCategory}},
  from a list of gene identifiers and their annotated GO categories.
  For each of the GO categories stated in \code{categ}, it also
  considers the implied annotation with all ancestors.
}
\usage{
cateGOryMatrix(x, categ)
}
\arguments{
  \item{x}{Character vector with (arbitrary) gene identifiers. They will
  be used for the column names of the resulting matrix.}
\item{categ}{A character vector of the same length as \code{x}
    with GO annotations for the genes in \code{x}. If a
    gene has multiple GO annotations, it is expected to occur multiple
    times in \code{x}, once for each different annotation.}
}
\details{
  Requires the \code{\link[GO:GO-package]{GO}} package.
}
\value{
 An integer matrix whose columns correspond to the unique elements of
 \code{x} and whose rows correspond to GO categories. Row names are the
 GO category names, column names the gene identifiers. Matrix element
 \code{[i,j]} is 1 if gene \code{j} is in category \code{i}, and 0 otherwise.
}

\author{W. Huber}

\seealso{  \code{\link{applyByCategory}} }
\examples{
  cateGOryMatrix(c("CG2671", "CG2671", "CG2950"),
                 c("GO:0000074", "GO:0001738", "GO:0003676"))
}

\keyword{manip}

