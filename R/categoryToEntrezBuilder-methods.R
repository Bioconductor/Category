setMethod("categoryToEntrezBuilder",
          signature(p="GeneKeggHyperGeoTestParams"),
          function(p) {
              keep.all <- switch(p@testDirection,
                                 over=FALSE,
                                 under=TRUE,
                                 stop("Bad testDirection slot"))
              getKeggToEntrezMap(p@geneIds, p@annotation, NULL, 
                                 p@universeGeneIds,
                                 keep.all=keep.all)
          })

setMethod("categoryToEntrezBuilder",
          signature(p="GeneGoHyperGeoTestParams"),
          function(p) {
              keep.all <- switch(p@testDirection,
                                 over=FALSE,
                                 under=TRUE,
                                 stop("Bad testDirection slot"))
              getGoToEntrezMap(p@geneIds, p@annotation, p@ontology, 
                               p@universeGeneIds, keep.all=keep.all)
          })


getGoToEntrezMap <- function(selected, lib, ontology, universe,
                             keep.all) {
    ## Return a list mapping GO ids to the Entrez Gene ids annotated
    ## at the GO id.  Only those GO ids that are in the specified
    ## ontology and have at least one annotation in the set of 
    ## Entrez Gene ids specified by 'selected' are included.
    go2allprobes <- getDataEnv("GO2ALLPROBES", lib)
    probeAnnot <- getGoToProbeMap(go2allprobes, ontology)
    ## Map to Entrez Gene and flag GO ids that don't have any
    ## annotations in our selected set.  No sense testing these.
    probeToEntrezMapHelper(probeAnnot, selected, lib, universe,
                           keep.all=keep.all)
}


getKeggToEntrezMap <- function(selected, lib, ontology, universe,
                               keep.all) {
    kegg2allprobes <- getDataEnv("PATH2PROBE", lib)
    probeAnnot <- getKeggToProbeMap(kegg2allprobes)
    probeToEntrezMapHelper(probeAnnot, selected, lib, universe,
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
    egAnnot <- lapply(probeAnnot, function(x) {
        ## YEAST doesn't use Entrez Gene, everybody else does
        z <- unique(x)
        if (lib != "YEAST")
          z <- unique(unlist(mget(unique(x), getDataEnv("LOCUSID", lib))))
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
    ## The elements are vectors of Probe probe ids.  Names are GO ids.
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
    inOnt = sapply(mget(goids, GOTERM), function(x) Ontology(x) == ontology)
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


removeLengthZeroAndMissing <- function(map) {
    notNA = sapply(map, function(x) {
        len = length(x)
        !(len == 0 || (len == 1 && is.na(x)))
    })
    map <- map[notNA]
}    

