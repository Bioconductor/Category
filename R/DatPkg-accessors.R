setMethod("ID2GO", "DatPkg",
          function(p) getAnnMap("GO", p@name))


setMethod("ID2GO", "OrganismMappingDatPkg",
          function(p) {
              ## FIXME: This is REALLY ugly
              ## The data in the *Mapping packages should
              ## be more sane so we don't have to munge it.
              e <- getAnnMap("LL2GO", chip=p@name)
              ans <- eapply(e, function(x) {
                  sub("([^@]+)@[A-Z]{3}", "\\1", x)
              })
              l2e(ans)
              })


setMethod("ID2EntrezID", "AffyDatPkg",
          function(p) getAnnMap("ENTREZID", p@name))

.createIdentifyMap <- function(keys) {
    e <- new.env(parent=emptyenv(), hash=TRUE)
    for (n in keys) {
        e[[n]] <- n
    }
    e
}

setMethod("ID2EntrezID", "YeastDatPkg",
          function(p) {
              .createIdentifyMap(ls(getAnnMap("CHR", p@name)))
          })

setMethod("ID2EntrezID", "Org.XX.egDatPkg",
          function(p) {
              .createIdentifyMap(ls(getAnnMap("CHR", p@name)))
          })


setMethod("ID2EntrezID", "OrganismMappingDatPkg",
          function(p) {
              ## Create an identity map
              ## FIXME: I'm not sure here if we really need the
              ## list of _all_ Entrez IDs or just those that have
              ## a GO mapping.  I think, in the end, that is all
              ## we need.  If this isn't enough, then we should
              ## take the unique(c(x1, x2, x2)) where the x_i
              ## are ls(humanLLMappingsLL2ACCNUM, 2GO, 2UG).
              e <- new.env(parent=emptyenv(), hash=TRUE)
              for (n in ls(getAnnMap("LL2GO", p@name))) {
                  e[[n]] <- n
              }
              e
          })


setMethod("GO2AllProbes", "DatPkg",
          function(p, ontology=c("BP", "CC", "MF")) {
              ontIds <- aqListGOIDs(ontology)
              go2all <- getAnnMap("GO2ALLPROBES", p@name)
              go2allOnt <- mget(ontIds, go2all, ifnotfound=NA)
              go2allOnt <- removeLengthZeroAndMissing(go2allOnt)
              l2e(go2allOnt)
          })


setMethod("GO2AllProbes", "OrganismMappingDatPkg",
          function(p, ontology=c("BP", "CC", "MF")) {
              ontIds <- aqListGOIDs(ontology)
              go2eg <- mget(ontIds, getAnnMap("GO2LL", p@name),
                            ifnotfound=NA)
              go2eg <- l2e(removeLengthZeroAndMissing(go2eg))
              goEnvName <- paste(ontology, "OFFSPRING", sep="")
              offspring <- mget(ls(go2eg), getAnnMap(goEnvName, "GO"),
                                ifnotfound=NA)
              go2allEg <- new.env(parent=emptyenv(), hash=TRUE,
                                  size=length(go2eg)*1.20)
              for (goid in names(offspring)) {
                  goids <- c(goid, offspring[[goid]])
                  goids <- goids[!is.na(goids)]
                  if (length(goids)) {
                      egids <- mget(goids, go2eg, ifnotfound=NA)
                      egids <- unique(unlist(egids))
                      go2allEg[[goid]] <- egids[!is.na(egids)]
                  }
              }
              go2allEg
          })

setMethod("GO2AllProbes", "Org.XX.egDatPkg",
          function(p, ontology=c("BP", "CC", "MF")) {

              db <- get("db_conn", paste("package:", p@name, sep=""))
              sqlQ <- "SELECT DISTINCT _left.gene_id,_right.go_id
              FROM genes AS _left INNER JOIN go_%s AS _right
              ON _left.id=_right.id WHERE _right.go_id IS NOT NULL"
              sqlQ <- sprintf(sqlQ, tolower(ontology))
              go2egTab <- dbGetQuery(db, sqlQ)
              go2eg <- l2e(split(go2egTab[["gene_id"]], go2egTab[["go_id"]]))

              goEnvName <- paste(ontology, "OFFSPRING", sep="")
              offspring <- mget(ls(go2eg),
                                getAnnMap(goEnvName, "GO"),
                                ifnotfound=NA)
              go2allEg <- new.env(parent=emptyenv(), hash=TRUE,
                                  size=length(go2eg)*1.20)
              for (goid in names(offspring)) {
                  goids <- c(goid, offspring[[goid]])
                  goids <- goids[!is.na(goids)]
                  if (length(goids)) {
                      egids <- mget(goids, go2eg, ifnotfound=NA)
                      egids <- unique(unlist(egids))
                      go2allEg[[goid]] <- egids[!is.na(egids)]
                  }
              }
              go2allEg
          })
