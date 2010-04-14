setMethod("universeBuilder", signature(p="KEGGHyperGParams"),
          function(p) {
                getUniverseViaKegg(p)
          })

setMethod("universeBuilder", signature(p="GOHyperGParams"),
          function(p) {
                getUniverseViaGo(p)
          })

setMethod("universeBuilder", signature(p="PFAMHyperGParams"),
          function(p) {
                getUniverseViaPfam(p)
          })

##The _db versions of these seem to never have been finished or been used.
XXXgetUniverseViaGo_db <- function(p) {
    datPkg <- p@datPkg
    ontology <- ontology(p)
    entrezIds <- universeGeneIds(p)
    ## Return all Entrez Gene Ids that are annotated at one or more
    ## GO terms belonging to the specified GO ontology.
    ## If 'entrezIds' is given, return the intersection of 'entrezIds'
    ## and the normal return value.
    ontology <- match.arg(ontology, c("BP", "CC", "MF"))
    ## FIXME: this might be faster as full-column pull and in-memory join
    SQL <- "select distinct gene_id from genes, go_XX where genes._id = go_XX._id"
    SQL <- gsub("XX", ontology, SQL)
    univ <- dbGetQuery(p@datPkg@getdb(), SQL)[[1]]
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
    ## FIXME: put dispatch code here depending on whether we get DB-based
    ## maps or env-based maps
    ontIds <- aqListGOIDs(ontology)
    probe2go <- eapply(ID2GO(datPkg), function(goids) {
        if (length(goids) == 0 || (length(goids) == 1 && is.na(goids)))
          return(FALSE)
        ## FIXME: *Mapping packages don't have a list,
        ##        but just the GO ID, so we have to branch
        if (!is.character(goids[[1]])) ## the Affy case
          goids <- subListExtract(goids, "GOID", simplify=TRUE)
        if (any(goids %in% ontIds))
          return(TRUE)
        FALSE
    })
    probe2go <- probe2go[unlist(probe2go)]
    probes <- names(probe2go)
    getUniverseHelper(probes, datPkg, entrezIds)
}


getUniverseViaKegg_db <- function(p) {
    entrezIds <- universeGeneIds(p)
    SQL <- "select distinct gene_id from genes, kegg where genes._id = kegg._id"
    univ <- dbGetQuery(p@datPkg@getdb(), SQL)[[1]]
    if (!is.null(entrezIds) && length(entrezIds) > 0)
      univ <- intersect(univ, unlist(entrezIds))
    if (length(univ) < 1)
      stop("No Entrez Gene ids left in universe")
    univ
}


getUniverseViaKegg <- function(p) {
    entrezIds <- universeGeneIds(p)
##    probe2kegg <- as.list(getDataEnv("PATH", annotation(p)))
    probe2kegg <- as.list(ID2KEGG(p@datPkg))
    notNA <- sapply(probe2kegg, function(x) !(length(x) == 1 && is.na(x)))
    probe2kegg <- probe2kegg[notNA]
    probes <- names(probe2kegg)
    getUniverseHelper(probes, p@datPkg, universeGeneIds(p))
}


getUniverseViaPfam_db <- function(p) {
    entrezIds <- universeGeneIds(p)
    SQL <- "select distinct gene_id from genes, pfam where genes._id = pfam._id"
    univ <- dbGetQuery(p@datPkg@getdb(), SQL)[[1]]
    if (!is.null(entrezIds) && length(entrezIds) > 0)
      univ <- intersect(univ, unlist(entrezIds))
    if (length(univ) < 1)
      stop("No ids left in gene universe")
    univ
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
    if (length(univ) < 1)
      stop("After filtering, there are no valid IDs that can be used as the Gene universe.\n  Check input values to confirm they are the same type as the central ID used by your annotation package.\n  For chip packages, this will still mean the central GENE identifier used by the package (NOT the probe IDs).")
    univ
}

