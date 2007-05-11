\name{HyperGResult-accessors}
\alias{HyperGResult-accessors}

\alias{pvalues}
\alias{pvalues,HyperGResult-method}
\alias{pvalues,ChrMapHyperGResult-method}

\alias{geneCounts}
\alias{geneCounts,HyperGResultBase-method}

\alias{universeCounts}
\alias{universeCounts,HyperGResultBase-method}

\alias{universeMappedCount}
\alias{universeMappedCount,HyperGResultBase-method}

\alias{geneMappedCount}

\alias{chrGraph}
\alias{chrGraph,ChrMapHyperGResult-method}

\alias{annotation}

\alias{description}
\alias{description,HyperGResultBase-method}

\alias{annotation,HyperGResultBase-method}

\alias{geneIds}
\alias{geneIds,HyperGResultBase-method}

\alias{geneIdUniverse}
\alias{geneIdUniverse,HyperGResult-method}
\alias{geneIdUniverse,ChrMapHyperGResult-method}

\alias{condGeneIdUniverse}
\alias{condGeneIdUniverse,HyperGResultBase-method}
\alias{condGeneIdUniverse,ChrMapHyperGResult-method}

\alias{geneIdsByCategory}
\alias{geneIdsByCategory,HyperGResultBase-method}

\alias{sigCategories}
\alias{sigCategories,HyperGResultBase-method}

\alias{geneMappedCount}
\alias{geneMappedCount,HyperGResultBase-method}

\alias{testName}
\alias{testName,HyperGResultBase-method}

\alias{isConditional,HyperGResultBase-method}
\alias{isConditional,ChrMapHyperGResult-method}

\alias{pvalueCutoff}
\alias{pvalueCutoff,HyperGResultBase-method}
\alias{testDirection}
\alias{testDirection,HyperGResultBase-method}


\alias{oddsRatios}
\alias{oddsRatios,HyperGResult-method}
\alias{oddsRatios,ChrMapHyperGResult-method}

\alias{expectedCounts}
\alias{expectedCounts,HyperGResult-method}
\alias{expectedCounts,ChrMapHyperGResult-method}

\alias{summary,HyperGResultBase-method}

\alias{htmlReport}
\alias{htmlReport,HyperGResultBase-method}


\docType{genericFunction}

\title{Accessors for HyperGResult Objects}

\description{
  This manual page documents generic functions for extracting data
  from the result object returned from a call to \code{hyperGTest}.
  The result object will be a subclass of \code{HyperGResultBase}.
  Methods apply to all result object classes unless otherwise noted.
}

\usage{

pvalues(r)
oddsRatios(r)
expectedCounts(r)

geneCounts(r)
universeCounts(r)
universeMappedCount(r)
geneMappedCount(r)

geneIds(r)
geneIdUniverse(r)
condGeneIdUniverse(r)
geneIdsByCategory(r, catids = NULL)
sigCategories(r, p)

## R CMD check doesn't like these
## annotation(r)
## description(r)

testName(r)
pvalueCutoff(r)
testDirection(r)

chrGraph(r)

}

\arguments{
  \item{r}{An instance of a subclass of
    \code{HyperGResultBase}.}

  \item{catids}{A character vector of category identifiers.}
  \item{p}{Numeric p-value used as a cutoff for selecting a subset of
    the result.}
}

\section{Accessor Methods (Generic Functions)}{
 \describe{

  \item{geneCounts}{returns an \code{"integer"}
      vector: for each category term tested, the number of genes from
      the gene set that are annotated at the term.}

  \item{pvalues}{returns a \code{"numeric"}
      vector: the ordered p-values for each category term tested.}

  \item{universeCounts}{returns an \code{"integer"}
      vector: for each category term tested, the number of genes from
      the gene universe that are annotated at the term.}

  \item{universeMappedCount}{returns an \code{"integer"}
      vector of length one giving the size of the gene universe set.}

  \item{expectedCounts}{returns a \code{"numeric"} vector
      giving the expected number of genes in the selected gene list to
      be found at each tested category term.}

  \item{oddsRatios}{returns a \code{"numeric"} vector giving
      the odds ratio for each category term tested.}

    \item{annotation}{returns the name
        of the annotation data package used.  }

    \item{geneIds}{returns the input
        vector of gene identifiers intersected with the universe of
        gene identifiers used in the computation. }

    \item{geneIdUniverse}{returns a list named by the
      tested categories.  Each element of the list is a vector of gene
      identifiers (from the gene universe) annotated at the
      corresponding category term.}

    \item{geneIdsByCategory}{returns a list similar to
      \code{geneIdUniverse}, but each vector of gene IDs is
      intersected with the list of selected gene IDs from
      \code{geneIds}.  The result is the selected gene IDs annotated at
      each category.}

    \item{sigCategories}{returns a character vector of category
      identifiers with a significant p-value.  If argument \code{p} is
      missing, then the cutoff obtained from \code{pvalueCutoff(r)} will
    be used.}

    \item{geneMappedCount}{returns the size of
        the selected gene set used in the computation.  This is simply
        \code{length(geneIds(obj))}.}

    \item{pvalueCutoff}{accessor for the
        \code{pvalueCutoff} slot.}

    \item{testDirection}{accessor for the
        \code{testDirection} slot.  Contains a string indicating
        whether the test was for \code{"over"} or \code{"under"}
        representation of the categories.}

    \item{description}{returns a character
        string description of the test result. }

    \item{testName}{returns a string
        describing the testing method used.}

    \item{isConditional}{returns \code{TRUE} if the
          result was obtained using a conditional algorithm.}

    \item{summary}{returns
      a \code{data.frame} summarizing the test result.  Optional
      arguments \code{pvalue} and \code{categorySize} allow
      specification of maximum p-value and minimum cateogyrSize,
      respectively.}

    \item{htmlReport}{writes
      an HTML version of the table produced by the
      \code{summary} method.  The path of a file to write the
      report to can be specified using the \code{file} argument.
      The default is \code{file=""} which will cause the report to
      be printed to the screen.  If you wish to create a single
      report comprising multiple results you can set
      \code{append=TRUE}.  The default is \code{FALSE} (overwrite
      pre-existing report file).  You can specify a string to use
      as an identifier for each table by providing a value for the
      \code{label} argument.  Additional named arguments will be
      passed to the \code{summary} method.}

  }
}

\details{
}

\author{Seth Falcon}
\seealso{
  \code{\link{hyperGTest}}
  \code{\link{HyperGResult-class}}
  \code{\link{HyperGParams-class}}
  \code{\link{GOHyperGParams-class}}
  \code{\link{KEGGHyperGParams-class}}
}
\keyword{htest}
