\name{ChrBandTree-class}
\docType{class}
\alias{ChrBandTree-class}

\alias{allGeneIds}
\alias{childrenOf}
\alias{parentOf}
\alias{lgeneIds}
\alias{treeLevels}
\alias{level2nodes}

\alias{allGeneIds,ChrBandTree-method}
\alias{childrenOf,ChrBandTree,character-method}
\alias{geneIds,ChrBandTree-method}
\alias{lgeneIds,ChrBandTree-method}
\alias{parentOf,ChrBandTree,character-method}
\alias{treeLevels,ChrBandTree-method}
\alias{level2nodes,ChrBandTree,numeric-method}
\alias{level2nodes,ChrBandTree,character-method}

\alias{show,ChrBandTree-method}

\title{Class "ChrBandTree"}
\description{
  This class represents chromosome band annotation data for a
given experiment.  The class is responsible for storing the mapping of
band to set of gene IDs located within that band as well as for
representing the tree structured relationship among the bands.

}

\section{Objects from the Class}{
Objects should be created using \code{NewChrBandTree} or
\code{ChrBandTreeFromGraph}.
}
\section{Slots}{
	 \describe{
    \item{\code{toParentGraph}:}{Object of class \code{"graph"}
    representing the tree of chromosome bands.  Edges in this directed
    graph go from child to parent.}
    \item{\code{toChildGraph}:}{Object of class \code{"graph"}.  This is
    the same as \code{toParentGraph}, but with the edge directons
    reversed.  This is not an ideal implementation due to the
    duplication of data, but it provides quick access to parents or
    children of a given node.}
    \item{\code{root}:}{Object of class \code{"character"} giving the
    name of the root node.  The convention is to use "ORGANISM:<organism>".}
    \item{\code{level2nodes}:}{Object of class \code{"list"} providing a
    mapping of levels in the tree to the set of nodes at that level.
    Levels \code{X} is defined as the set of nodes with a path length of
    \code{X} from the root node.}
  }
}
\section{Methods}{
  \describe{

    \item{allGeneIds}{Return a vector of gene IDs representing the gene
    universe for this \code{ChrBandTree}}

  \item{childrenOf}{Return a list with an element for each the
    character vector \code{n}.  Each element is a character vector of
    node names of the children of the named element.}

  \item{geneIds}{Return a vector of gene IDs for a single band.}
  
  \item{lgeneIds}{Return a list of vectors of gene IDs when given more
    than one band.  The "l" prefix is for list.}
  
  \item{parentOf}{Return the parents of the specified bands.  See
    \code{childrenOf} for a description of the structure of the return
    value.}

  \item{treeLevels}{Return an integer vector identifying the levels of
    the tree.}

  \item{level2nodes(g, level)}{Return the nodes in the tree that are at
    the level specified by \code{level}.  The \code{level} argument can
    be either numeric or character, but should match a level returned by
    \code{treeLevels}.}
}
}
\author{S. Falcon}
\note{
  Not all known chromosome bands will be represented in a given
  instance.  The set of bands that will be present is determined by the
  available annotation data and the specified gene universe.  The
  annotation source maps genes to their most specific band.  Such bands
  and all bands on the path to the root will be represented in the
  resulting tree.

  Currently there is only support for human and mouse data.
}

\examples{
library("hgu95av2.db")
set.seed(0xfeee)
univ = NULL ## use all Entrez Gene IDs on the chip (not recommended)
ct = NewChrBandTree("hgu95av2.db", univ)

length(allGeneIds(ct))

exampleLevels(ct)

geneIds(ct, "10p11")
lgeneIds(ct, "10p11")
lgeneIds(ct, c("10p11", "Yq11.22"))

pp = parentOf(ct, c("10p11", "Yq11.22"))
childrenOf(ct, unlist(pp))

treeLevels(ct)

level2nodes(ct, 0)
level2nodes(ct, 0L)
level2nodes(ct, "0")

level2nodes(ct, 1)

}
\keyword{classes}
