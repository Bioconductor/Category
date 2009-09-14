\name{cateGOry}
\alias{cateGOry}
\title{Construct a category membership matrix from a list of gene
  identifiers and their annotated GO categories.}
\description{
  The function constructs a category membership matrix, such as used by
  \code{\link{applyByCategory}},
  from a list of gene identifiers and their annotated GO categories.
  For each of the GO categories stated in \code{categ},
  all less specific terms (ancestors) are also included, thus one need
  only obtain the most specific set of GO term mappings, which
  can be obtained from Bioconductor annotation packages or via \pkg{biomaRt}.
  The ancestor relationships are obtained from the \pkg{GO.db} package.   
}
\usage{
cateGOry(x, categ, sparse=FALSE)
}
\arguments{
  \item{x}{Character vector with (arbitrary) gene identifiers. They will
  be used for the column names of the resulting matrix.}
\item{categ}{A character vector of the same length as \code{x}
    with GO annotations for the genes in \code{x}. If a
    gene has multiple GO annotations, it is expected to occur multiple
    times in \code{x}, once for each different annotation.}
  \item{sparse}{Logical. Currently, this is ignored. 
    This argument might be used in future versions of the function to
    result in returning a sparse matrix representation.}
}
\details{
  Requires the \code{\link[GO.db:GOBASE]{GO}} package.

  For some subsequent analyses, it is useful to remove categories that
  have only a small number of members. Use the normal matrix subsetting
  syntax for this, see example.

  If a GO category in \code{categ} is not found in the GO annotation
  package, a warning will be generated, and no ancestors
  for that GO category are added (but that category itself will be part
  of the returned adjacency matrix).
}
\value{
  The adjacency matrix of the bipartite category membership graph,
  rows are categories and columns genes.
}

\author{W. Huber}

\seealso{\code{\link{applyByCategory}}}
\examples{
  g = cateGOry(c("CG2671", "CG2671", "CG2950"),
               c("GO:0000074", "GO:0001738", "GO:0003676"))
  g

  rownames(g)
  colnames(g)
  rowSums(g)   ## number of genes in each category

  ## Filter out categories with less than minMem members.
  ## This is toy data, in real applications, a number higher
  ## than 2 will be more appropriate.
  minMemb = 2
  g[rowSums(g)>=minMemb,,drop=FALSE ]
}

\keyword{manip}

