setMethod("categoryToEntrezBuilder",
          signature(p="KEGGHyperGParams"),
          function(p) {
              getKeggToEntrezMap(p)
          })

setMethod("categoryToEntrezBuilder",
          signature(p="GOHyperGParams"),
          function(p) {
              getGoToEntrezMap(p)
          })

setMethod("categoryToEntrezBuilder",
          signature(p="PFAMHyperGParams"),
          function(p) {
              getPfamToEntrezMap(p)
          })


getGoToEntrezMap <- function(p) {
    keep.all <- switch(testDirection(p),
                       over=FALSE,
                       under=TRUE,
                       stop("Bad testDirection slot"))    
    lib <- p@datPkg
    ontology <- ontology(p)
    ## Return a list mapping GO ids to the Entrez Gene ids annotated
    ## at the GO id.  Only those GO ids that are in the specified
    ## ontology and have at least one annotation in the set of 
    ## Entrez Gene ids specified by 'selected' are included.
    go2allprobes <- GO2AllProbes(lib, ontology)
    probeAnnot <- getGoToProbeMap(go2allprobes, ontology)
    ## Map to Entrez Gene and flag GO ids that don't have any
    ## annotations in our selected set.  No sense testing these.
    probeToEntrezMapHelper(probeAnnot, geneIds(p), lib, universeGeneIds(p),
                           keep.all=keep.all)
}


getKeggToEntrezMap <- function(p) {
    keep.all <- switch(testDirection(p),
                       over=FALSE,
                       under=TRUE,
                       stop("Bad testDirection slot"))
    lib <- annotation(p)
    kegg2allprobes <- getDataEnv("PATH2PROBE", lib)
    probeAnnot <- getKeggToProbeMap(kegg2allprobes)
    probeToEntrezMapHelper(probeAnnot, geneIds(p), p@datPkg, universeGeneIds(p),
                           keep.all=keep.all)
}


getPfamToEntrezMap <- function(p) {
    keep.all <- switch(testDirection(p),
                       over=FALSE,
                       under=TRUE,
                       stop("Bad testDirection slot"))
    probe2pfam <- getDataEnv("PFAM", annotation(p))
    pfam2allprobes <- splitOrfByPfam(probe2pfam)
    probeAnnot <- getPfamToProbeMap(pfam2allprobes)
    probeToEntrezMapHelper(probeAnnot, geneIds(p), p@datPkg, universeGeneIds(p),
                           keep.all=keep.all)
}


probeToEntrezMapHelper <- function(probeAnnot, selected, lib, universe,
                                   keep.all=FALSE) {
    ## Given a list 'probeAnnot' mapping category => probe, convert the probe
    ## ids to Entrez Gene ids (unless we are using YEAST, in which case we skip
    ## this step).  Then reduce the Entrez Gene ids to the specified universe
    ## and only keep those entries in the list that have a non-empty
    ## intersection with the selected genes.  If keep.all is TRUE, then we keep
    ## entries even if the list of gene IDs includes no gene ID from the
    ## selected list.
    id2entrezEnv <- ID2EntrezID(lib)
    egAnnot <- lapply(probeAnnot, function(x) {
        z <- unique(x)
        z <- unique(unlist(mget(unique(x), id2entrezEnv)))
        z  <- intersect(z, universe)
        ## would be nice to have a short-circuiting way to do this
        if (length(z) > 0 && (keep.all || any(selected %in% z))) {
              return(z)
        }
        NULL
    })
    notNull <- sapply(egAnnot, function(x) !is.null(x))
    egAnnot[notNull]
}


getGoToProbeMap <- function(go2allprobes, ontology, goids) {
    ## Return a list with one element for each GO id in the specified 
    ## ontology that has at least one Probe probe id annotated at it.  
    ## The elements are vectors of Probe ids.  Names are GO ids.
    ##
    probeAnnot = as.list(go2allprobes)
    if (!missing(goids))
      probeAnnot = probeAnnot[goids]

    ## Remove "all" artifact GO nodes
    whAll = match("all", names(probeAnnot), 0)
    ## This is a tad faster than if (any(whAll > 0)) 
    ## because we stop early
    for (el in whAll) {
        if (el > 0) {
            probeAnnot = probeAnnot[-whAll]
            break
        }
    }

    goids = names(probeAnnot)

    ## filter on desired GO ontology
    inOnt <- filterGOByOntology(goids, ontology)
    probeAnnot = probeAnnot[inOnt]

    ## remove any GO ids that don't map to an probe id (NA)
    removeLengthZeroAndMissing(probeAnnot)
}


getKeggToProbeMap <- function(kegg2allprobes, keggIds) {
    probeAnnot = as.list(kegg2allprobes)
    if (!missing(keggIds))
      probeAnnot = probeAnnot[keggIds]
    removeLengthZeroAndMissing(probeAnnot)
}

getPfamToProbeMap <- function(pfam2allprobes, pfamIds) {
    probeAnnot = as.list(pfam2allprobes)
    if (!missing(pfamIds))
      probeAnnot = probeAnnot[pfamIds]
    removeLengthZeroAndMissing(probeAnnot)
}



removeLengthZeroAndMissing <- function(map) {
    notNA = sapply(map, function(x) {
        len = length(x)
        !(len == 0 || (len == 1 && is.na(x)))
    })
    map <- map[notNA]
}    

splitOrfByPfam <- function(ypfEnv){
  ## There should be a PFAM to ORF data set in addition to YEASTPFAM...a
  ## YEASTPFAM2PROBE, but for now we invert
  probe2pfam <- as.list(ypfEnv)
  probe2pfam <- removeLengthZeroAndMissing(probe2pfam)
  pfamlen <- listLen(probe2pfam)
  orf <- rep(names(probe2pfam), pfamlen)
  pfam <- unlist(probe2pfam)
  orfByPfam <- split(orf, pfam)
  orfByPfam

}
