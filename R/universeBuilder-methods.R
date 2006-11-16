setMethod("universeBuilder", signature(p="KEGGHyperGParams"),
          function(p) {
              getUniverseViaKegg(p)
          })

setMethod("universeBuilder", signature(p="GOHyperGParams"),
          function(p) {
              ## FIXME: this should be done via dispatch
              if (is(p@datPkg, "DBPkg"))
                getUniverseViaGo_db(p)
              else
                getUniverseViaGo(p)
          })

setMethod("universeBuilder", signature(p="PFAMHyperGParams"),
          function(p) {
            getUniverseViaPfam(p)
          })


setGeneric("getUniverse", function(datPkg) standardGeneric("getUniverse"))


getUniverseViaGo_db <- function(p) {
    datPkg <- p@datPkg
    ontology <- ontology(p)
    entrezIds <- universeGeneIds(p)
    ## Return all Entrez Gene Ids that are annotated at one or more
    ## GO terms belonging to the specified GO ontology.
    ## If 'entrezIds' is given, return the intersection of 'entrezIds'
    ## and the normal return value.
    ontology <- match.arg(ontology, c("BP", "CC", "MF"))
    SQL <- "select distinct gene_id from GO, probe2gene where GO.Ontology ='%s' and GO.PROBE = probe2gene.probe_id"
    univ <- dbGetQuery(p@datPkg@getdb(), sprintf(SQL, ontology)[[1]]
    if (!is.null(entrezIds) && length(entrezIds) > 0)
      univ <- intersect(univ, unlist(entrezIds))
    if (length(univ) < 1)
      stop("No Entrez Gene ids left in universe")
    univ
}


getUniverseViaGo <- function(p) {
    datPkg <- p@datPkg
    ontology <- ontology(p)
    entrezIds <- universeGeneIds(p)
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

getUniverseViaKegg <- function(p) {
    entrezIds <- universeGeneIds(p)
    probe2kegg <- as.list(getDataEnv("PATH", annotation(p)))
    notNA <- sapply(probe2kegg, function(x) !(length(x) == 1 && is.na(x)))
    probe2kegg <- probe2kegg[notNA]
    probes <- names(probe2kegg)
    getUniverseHelper(probes, p@datPkg, universeGeneIds(p))
}


getUniverseViaPfam <- function(p) {
    entrezIds <- universeGeneIds(p)
    probe2pfam <- as.list(getDataEnv("PFAM", annotation(p)))
    notNA <- sapply(probe2pfam, function(x) !(length(x) == 1 && is.na(x)))
    probe2pfam <- probe2pfam[notNA]
    probes <- names(probe2pfam)
    getUniverseHelper(probes, p@datPkg, entrezIds)
}



getUniverseHelper <- function(probes, datPkg, entrezIds) {
    univ <- unique(unlist(mget(probes, ID2EntrezID(datPkg))))
    if (!missing(entrezIds) && !is.null(entrezIds) && length(entrezIds) > 0)
      univ <- intersect(univ, unlist(entrezIds))
    if (length(univ) < 1) ##FIXME: improve error msg
      stop("No Entrez Gene ids left in universe")
    univ
}

