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


getGoToEntrezMap_db <- function(p, db) {
    keep.all <- switch(testDirection(p), over=FALSE, under=TRUE,
                       stop("Bad testDirection slot"))    
    ##db <- p@datPkg@getdb()
    ontology <- ontology(p)
    GOXXALL <- sprintf("GO%sALL", ontology)
    goxxall.id <- sqliteQuickColumn(db, GOXXALL, "ID")
    goxxall.goid <- sqliteQuickColumn(db, GOXXALL, "GOID")
    genes.id <- sqliteQuickColumn(hgdb, "GENES", "ID")
    genes.geneid <- sqliteQuickColumn(db, "GENES", "GeneID")
    ## convert internal gene id to "official" gene id (usually Entrez)
    intId2Entrez <- split(genes.geneid, genes.id)
    ## FIXME: not all IDs will be integers!
    goxxall.geneid <- as.integer(unlist(intId2Entrez[goxxall.id]))
    go2all <- split(goxxall.geneid, goxxall.goid)
    sql <- paste("select distinct GOID from", GOXXALL,
                 ", GENES where GeneID in (",
                 paste("'", geneIds(p), "'", sep="", collapse=","), ")",
                 "and ", paste(GOXXALL, ".ID", sep=""), "= GENES.ID")
    rs <- dbSendQuery(db, sql)
    on.exit(dbClearResult(rs))
    ## FIXME: can we estimate size?
    goids <- fetch(rs, n=-1)[[1]]
    ## FIXME, why don't we get back a unique result set?
    goids <- unique(goids)
    stopifnot(dbHasCompleted(rs))
    univ <- universeGeneIds(p)
    go2all <- lapply(go2all[goids], function(eg) {
        z <- intersect(univ, eg)
        stopifnot(length(z) > 0)
        z
    })
    go2all
}
    

getGoToEntrezMap_db2 <- function(p, db) {
    keep.all <- switch(testDirection(p), over=FALSE, under=TRUE,
                       stop("Bad testDirection slot"))    
    ##db <- p@datPkg@getdb()
    ontology <- ontology(p)
    GOXXALL <- sprintf("GO%sALL", ontology)
    tables <- paste(GOXXALL, "GENES", sep=", ")
    ourGeneIds <- paste("'", geneIds(p), "'", sep="", collapse=",")
    sql <- paste("select distinct GOID from", tables,
                 "where GeneID in (", ourGeneIds, ")",
                 "and ", paste(GOXXALL, ".ID", sep=""), "= GENES.ID")
    sql <- paste("select GeneID, GOID from", tables,
                 "where GOID in (", sql, ")",
                 "and ", paste(GOXXALL, ".ID", sep=""), "= GENES.ID")
    rs <- dbSendQuery(db, sql)
    on.exit(dbClearResult(rs))
    ans <- fetch(rs, n=-1)
    go2all <- split(ans[["GeneID"]], ans[["GOID"]])
    ##SQL <- "select GeneID, GOBPALL.GOID from GOXXALL, GENES where GOXXALL.ID = GENES.ID"
    ##SQL <- gsub("XX", ontology, SQL)
    univ <- unlist(universeGeneIds(p), use.names=FALSE)
    go2all <- lapply(go2all, function(eg) {
        z <- intersect(univ, eg)
        stopifnot(length(z) > 0)
        z
    })
    go2all
}

    


##     sql <- paste("select distinct ID from GENES where GeneID in (",
##                  paste("'", myEg, "'", sep="", collapse=","), ")")
##     rs <- dbSendQuery(db, sql)
##     ans <- fetch(rs, n=length(myEg)+2L)[[1]]
##     stopifnot(dbHasCompleted(rs))
##     dbClearResult(rs)


##     ## This is slower than pasting in the IN (x, y, z, ...) clause
##     ST({
##         dbBeginTransaction(db)
##         sql1 <- "create temp table idinput (id int)"
##         dbGetQuery(db, sql1)
##         sql2 <- "insert into idinput values (?)"
##         dbGetPreparedQuery(db, sql2, data.frame(id=myEg))
##         sql3 <- paste("select distinct ID from GENES where GeneID in (",
##                      "select id from idinput)")
##         rs <- dbSendQuery(db, sql3)
##         ans <- fetch(rs, n=length(myEg)+2L)
##         stopifnot(dbHasCompleted(rs))
##         dbClearResult(rs)
##         dbGetQuery(db, "drop table idinput")
##         dbCommit(db)
##         ans <- ans[[1]]
##     })
    



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
