setMethod("universeBuilder", signature(p="GeneKeggHyperGeoTestParams"),
          function(p) {
            getUniverseViaKegg(p@annotation, p@universeGeneIds)  
          })

setMethod("universeBuilder", signature(p="GeneGoHyperGeoTestParams"),
          function(p) {
            getUniverseViaGo(p@annotation, p@ontology, p@universeGeneIds)
        })


getUniverseViaGo <- function(lib, ontology="BP", entrezIds=NULL) {
    ## Return all Entrez Gene Ids that are annotated at one or more
    ## GO terms belonging to the specified GO ontology.
    ## If 'entrezIds' is given, return the intersection of 'entrezIds'
    ## and the normal return value.
    ontology <- match.arg(ontology, c("BP", "CC", "MF"))
    probe2go <- eapply(getDataEnv("GO", lib), function(goids) {
        if (length(goids) == 0 || (length(goids) == 1 && is.na(goids)))
          return(FALSE)
        ## Normally, would do sapply here, but we only want to know if
        ## at least one is TRUE.  It is a lot faster to stop short.
        for (goid in goids) {
            ## use identical in case the GO data is borked
            ## and Ontology is NA
            if (identical(goid$Ontology, ontology))
                return(TRUE)
        }
        FALSE
    })
    probe2go <- probe2go[unlist(probe2go)]
    probes <- names(probe2go)
    getUniverseHelper(probes, lib, entrezIds)
}


getUniverseViaKegg <- function(lib, entrezIds, ...) {
    probe2kegg <- as.list(getDataEnv("PATH", lib))
    notNA <- sapply(probe2kegg, function(x) !(length(x) == 1 && is.na(x)))
    probe2kegg <- probe2kegg[notNA]
    probes <- names(probe2kegg)
    getUniverseHelper(probes, lib, entrezIds)
}


getUniverseHelper <- function(probes, lib, entrezIds) {
    if (lib == "YEAST")
      univ <- probes
    else
      univ <- unique(unlist(mget(probes, getDataEnv("LOCUSID", lib))))
    if (!missing(entrezIds) && !is.null(entrezIds) && length(entrezIds) > 0)
      univ <- intersect(univ, entrezIds)
    if (length(univ) < 1) ##FIXME: improve error msg
      stop("No Entrez Gene ids left in universe")
    univ
}

