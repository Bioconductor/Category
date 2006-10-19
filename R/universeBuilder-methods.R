setMethod("universeBuilder", signature(p="KEGGHyperGParams"),
          function(p) {
            getUniverseViaKegg(p@annotation, p@universeGeneIds)  
          })

setMethod("universeBuilder", signature(p="GOHyperGParams"),
          function(p) {
            getUniverseViaGo(p@datPkg, p@ontology, p@universeGeneIds)
        })

setMethod("universeBuilder", signature(p="PFAMHyperGParams"),
          function(p) {
            getUniverseViaPfam(p@annotation, p@universeGeneIds)
          })


getUniverseViaGo <- function(datPkg, ontology="BP", entrezIds=NULL) {
    ## Return all Entrez Gene Ids that are annotated at one or more
    ## GO terms belonging to the specified GO ontology.
    ## If 'entrezIds' is given, return the intersection of 'entrezIds'
    ## and the normal return value.
    ontology <- match.arg(ontology, c("BP", "CC", "MF"))
    ontIds <- getGOOntologyIDs(ontology)
    probe2go <- eapply(ID2GO(datPkg), function(goids) {
        if (length(goids) == 0 || (length(goids) == 1 && is.na(goids)))
          return(FALSE)
        ## FIXME: *Mapping packages don't have a list,
        ##        but just the GO ID, so we have to branch
        if (!is.character(goids[[1]])) ## the Affy case
          goids <- sapply(goids, function(x) x$GOID)
        if (any(goids %in% ontIds))
          return(TRUE)
        FALSE
    })
    probe2go <- probe2go[unlist(probe2go)]
    probes <- names(probe2go)
    getUniverseHelper(probes, datPkg, entrezIds)
}


getUniverseViaKegg <- function(lib, entrezIds) {
    probe2kegg <- as.list(getDataEnv("PATH", lib))
    notNA <- sapply(probe2kegg, function(x) !(length(x) == 1 && is.na(x)))
    probe2kegg <- probe2kegg[notNA]
    probes <- names(probe2kegg)
    getUniverseHelper(probes, lib, entrezIds)
}

getUniverseViaPfam <- function(lib, entrezIds) {
    probe2pfam <- as.list(getDataEnv("PFAM", lib))
    notNA <- sapply(probe2pfam, function(x) !(length(x) == 1 && is.na(x)))
    probe2pfam <- probe2pfam[notNA]
    probes <- names(probe2pfam)
    getUniverseHelper(probes, lib, entrezIds)
}



getUniverseHelper <- function(probes, datPkg, entrezIds) {
    univ <- unique(unlist(mget(probes, ID2EntrezID(datPkg))))
    if (!missing(entrezIds) && !is.null(entrezIds) && length(entrezIds) > 0)
      univ <- intersect(univ, unlist(entrezIds))
    if (length(univ) < 1) ##FIXME: improve error msg
      stop("No Entrez Gene ids left in universe")
    univ
}

