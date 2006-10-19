setMethod("ID2GO", "DatPkg",
          function(p) getDataEnv("GO", p@name))


setMethod("ID2GO", "OrganismMappingDatPkg",
          function(p) getDataEnv("LL2GO", p@name))


setMethod("ID2EntrezID", "AffyDatPkg",
          function(p) getDataEnv("ENTREZID", p@name))

setMethod("ID2EntrezID", "YeastDatPkg",
          function(p) {
              ## Create an identity map
              e <- new.env(parent=emptyenv(), hash=TRUE)
              for (n in ls(getDataEnv("CHR", p@name))) {
                  e[[n]] <- n
              }
              e
          })


setMethod("ID2EntrezID", "OrganismMappingDatPkg",
          function(p) {
              ## Create an identity map
              e <- new.env(parent=emptyenv(), hash=TRUE)
              for (n in ls(getDataEnv("LL2ACCNUM", p@name))) {
                  e[[n]] <- n
              }
              e
          })


setMethod("GO2AllProbes", "DatPkg",
          function(p, ontology=c("BP", "CC", "MF")) {
              ontology <- match.arg(ontology)
              goEnvName <- paste(ontology, "PARENTS", sep="")
              ontIds <- ls(getDataEnv(goEnvName, "GO"))
              go2all <- getDataEnv("GO2ALLPROBES", p@name)
              go2allOnt <- mget(ontIds, go2all, ifnotfound=NA)
              go2allOnt <- removeLengthZeroAndMissing(go2allOnt)
              l2e(go2allOnt)
          })


setMethod("GO2AllProbes", "OrganismMappingDatPkg",
          function(p, ontology=c("BP", "CC", "MF")) {
              ontology <- match.arg(ontology)
              goEnvName <- paste(ontology, "PARENTS", sep="")
              ontIds <- ls(getDataEnv(goEnvName, "GO"))
              go2eg <- mget(ontIds, getDataEnv("GO2LL", p@name),
                            ifnotfound=NA)
              go2eg <- l2e(removeLengthZeroAndMissing(go2eg))
              goEnvName <- paste(ontology, "OFFSPRING", sep="")
              offspring <- mget(ls(go2eg), getDataEnv(goEnvName, "GO"),
                                ifnotfound=NA)
              go2allEg <- new.env(parent=emptyenv(), hash=TRUE)
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
              


##FIXME: add a getGOBPIDs, getGOCCIDs, getGOMFIDs to annotate
