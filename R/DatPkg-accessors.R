setMethod("ID2GO", "DatPkg",
          function(p) getAnnMap("GO", p@name))




setMethod("ID2EntrezID", "AffyDatPkg",
          function(p) getAnnMap("ENTREZID", p@name))

##FIXME: this is seriously slow - try l2e to speed up a bit
.createIdentityMap <- function(keys) {
    keys = as.list(keys)
    names(keys) = keys 
    l2e(keys)
#    e <- new.env(parent=emptyenv(), hash=TRUE)
#    for (n in keys) {
#        e[[n]] <- n
#    }
#    e
}

##this needs to handle all new, old and org based yeast packages
setMethod("ID2EntrezID", "YeastDatPkg",
          function(p) {
              bname = sub("\\.db$", "", p@name)
              if( exists( paste(bname, "ORF", sep="")) ) 
	        return(getAnnMap("ORF", p@name))
              else
              .createIdentityMap(ls(getAnnMap("CHR", p@name)))
          })

setMethod("ID2EntrezID", "ArabidopsisDatPkg",
          function(p) {
              bname = sub("\\.db$", "", p@name)
              if( exists( paste(bname, "ACCNUM", sep="")) ) 
	        return(getAnnMap("ACCNUM", p@name))
              else
              .createIdentityMap(ls(getAnnMap("CHR", p@name)))
          })

setMethod("ID2EntrezID", "Org.XX.egDatPkg",
          function(p) {
              .createIdentityMap(ls(getAnnMap("CHR", p@name)))
          })



setMethod("GO2AllProbes", "DatPkg",
          function(p, ontology=c("BP", "CC", "MF")) {
              ontIds <- aqListGOIDs(ontology)
              go2all <- getAnnMap("GO2ALLPROBES", p@name)
              ontIds <- intersect(ontIds, ls(go2all))
              go2allOnt <- mget(ontIds, go2all, ifnotfound=NA)
              go2allOnt <- removeLengthZeroAndMissing(go2allOnt)
              l2e(go2allOnt)
          })


setMethod("GO2AllProbes", "YeastDatPkg",
          function(p, ontology=c("BP", "CC", "MF")) {
              conn <- do.call(paste(sub("\\.db", "", p@name), "_dbconn", sep=""), list())
              schema <- dbmeta(conn, "DBSCHEMA")
              env = environment()
              if(schema == "YEASTCHIP_DB"){
                  env = callNextMethod()
                  return(env)
              }
              ontIds <- aqListGOIDs(ontology)
              go2all <- getAnnMap("GO2ALLORFS", p@name)
              ontIds <- intersect(ontIds, ls(go2all))
              go2allOnt <- mget(ontIds, go2all, ifnotfound=NA)
              go2allOnt <- removeLengthZeroAndMissing(go2allOnt)
              env = l2e(go2allOnt)
              return(env)
          })



setMethod("GO2AllProbes", "Org.XX.egDatPkg",
          function(p, ontology=c("BP", "CC", "MF")) {

              #db <- get("db_conn", paste("package:", p@name, sep=""))
              pname = sub("\\.db", "", p@name)
              db <- do.call(paste(pname, "dbconn", sep="_"), list())
              sqlQ <- "SELECT DISTINCT gene_id, go_id
              FROM genes INNER JOIN go_%s USING (_id)"
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

isDBDatPkg <- function(dpkg) {
    length(grep("\\.db$", dpkg@name)) > 0
}
