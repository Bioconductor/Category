\name{cb_contingency}
\alias{cb_contingency}
\alias{cb_sigBands}
\alias{cb_children}

\title{Create and Test Contingency Tables of Chromosome Band Annotations}
\description{
  For each chromosome band identifier in \code{chrVect},
  \code{cb_contingency} builds and performs a test on a 2 x k
  contingency table for the genes from \code{selids} found in the child
  bands of the given \code{chrVect} element.

  \code{cb_sigBands} extracts the chromosome band identifiers that were
  in a contingency table that tested significant given the specified
  p-value cutoff.

  \code{cb_children} returns the child bands of a given band in the
  chromosome band graph.  The argument must have length equal to one.
}
\usage{
cb_contingency(selids, chrVect, chrGraph, testFun = chisq.test,
               min.expected = 5L, min.k = 1L)

cb_sigBands(b, p.value = 0.01)

cb_children(n, chrGraph)
}

\arguments{
  \item{selids}{A vector of the selected gene identifiers (usual Entrez IDs).}
  \item{chrVect}{A character vector of chromosome band identifiers}
  \item{chrGraph}{A \code{graph} object as returned by
    \code{makeChrBandGraph}.  The nodes should be chromosome band IDs
    and the edges should represent the tree structure of the bands.
    Furthermore, the graph is expected to have a \code{"geneIds"} node
    attribute providing a vector of gene IDs annotated at each band.}
  \item{testFun}{The function to use for testing the 2 x k contingency
  tables.  The default is \code{chisq.test}.  It will be called with a
  single argument, a 2 x k matrix representing the contingency table.}
  \item{min.expected}{A numeric value specifying the minimum expected
    count for columns to be included in the contingency table.  The
    expected count is \code{(rowSum * colSum) / n}.  Chromosome bands
    with a select cell count less than \code{min.expected} are dropped
    from the table before testing occurs.}
  \item{min.k}{An integer giving the minimum number of chromosome bands
    that must be present in a contingency table in order to proceed with
    testing.}
  \item{b}{A list as returned by \code{cb_contingency}}
  \item{p.value}{A p-value cutoff to use in selecting significant
    contingency tables.}
  \item{n}{A length one character vector specifying a chromosome band
    annotation.  Bands not found in \code{chrGraph} will return
    \code{character(0)} when passed to \code{cb_children}.}

}

\details{
  \code{cb_sigBands} assumes that the p-value associated with a result
  of \code{testFun} can by accessed as \code{testFun(t)$p.value}.  We
  should improve this to be a method call which can then be specialized
  based on the class of the object returned by \code{testFun}.
}

\value{

  \code{cb_contingency} returns a list with an element for each test
  performed.  This will most often be shorter than
  \code{length(chrVect)} due to skipped tests based on \code{min.found}
  and \code{min.k}.  Each element of the returned list is itself a list
  with components:
  \item{table}{A 2 x k contingency table}
  \item{result }{The output of \code{testFun} applied to the \code{table}.}

  \code{cb_sigBands} returns a character vector of chromosome band
  identifiers that are in one of the contingency tables that had a
  p-value less than the cutoff specified by \code{p.value}.
}

\author{Seth Falcon}

\keyword{ }

