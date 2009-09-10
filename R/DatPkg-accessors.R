setMethod("ID2GO", "DatPkg",
          function(p) getAnnMap("GO", p@name))


setMethod("ID2GO", "GeneSetCollectionDatPkg",  ##yes, it is already a list, but its backwards of what we want here.
          function(p){
            coll <- p@GeneSetCollection
            genes <- geneIds(coll)
            genesLengths <- lapply(genes, length)
            GOIDs <- names(coll)
            GOIDReps <- rep(GOIDs, genesLengths)
            collFrame <- cbind(GOIDReps, unlist(genes))
            collList <- split(as.character(collFrame[,1]), as.character(collFrame[,2]))
            result <- l2e(collList)
            result
          })


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
              bname = p@name 
              if( exists( paste(bname, "ORF", sep="")) ) 
	        return(getAnnMap("ORF", p@name))
              else
              .createIdentityMap(allValidKeys(p@name))
          })

setMethod("ID2EntrezID", "ArabidopsisDatPkg",
          function(p) {
              bname = p@name 
              if( exists( paste(bname, "ACCNUM", sep="")) ) 
	        return(getAnnMap("ACCNUM", p@name))
              else
              .createIdentityMap(allValidKeys(p@name))
          })

setMethod("ID2EntrezID", "Org.XX.egDatPkg",
          function(p) {
              .createIdentityMap(allValidKeys(p@name))
          })



setMethod("ID2EntrezID", "GeneSetCollectionDatPkg",
          function(p) {
            ##This method does not need to really do anything "real" since
            ##they are going to get out the ID type that they put in: no
            ##matter what.
            coll <- p@GeneSetCollection
            genes <- unique(unlist(geneIds(coll)))
            collList <- split(genes,genes)
            res <- l2e(collList)
            res
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
              conn <- do.call(paste(p@name, "_dbconn", sep=""), list())
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
              db <- do.call(paste(p@name, "dbconn", sep="_"), list())
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



setMethod("GO2AllProbes", "GeneSetCollectionDatPkg",  
          function(p, ontology=c("BP", "CC", "MF")) {
            coll <- p@GeneSetCollection
            ## Lets put the GeneSetCollection into a format that is easier to
            ## filter (from the left OR right)
            genes = geneIds(coll)
            genesLengths = lapply(genes, length)
            GOIDs = names(coll)
            GOIDReps = rep(GOIDs, genesLengths)
            collFrame = cbind(GOIDReps, unlist(genes))
            
            ##Now filter out all GOIDs not from the selected ontology
            ontology <- ontology
            ontology <- match.arg(ontology, c("BP", "CC", "MF"))
            ontIds <- aqListGOIDs(ontology)
            ontFilt <- collFrame[,1] %in% ontIds
            collFrame <- collFrame[ontFilt,]
            
            ##Then put things back into a list format
            result <- split(as.character(collFrame[,2]), as.character(collFrame[,1]))
            
            if(length(result)==0){
              stop("Sorry, but there just aren't any annotations in your Gene Set Collection for any of the genes you want to test on.")
            }
            result <- l2e(result)
            result
          })



setMethod("isDBDatPkg","DatPkg",
          function(d){
            ##If there is a connection object then it's a db package.
            require(paste(d@name, ".db", sep=""), character.only=TRUE)
            exists(paste(d@name, "_dbconn", sep=""), mode="function")
          })


setMethod("isDBDatPkg","GeneSetCollectionDatPkg", function(d){return(FALSE)})

















####################################################################
## DatPkgFactory methods (default has to be back in AllClasses) :(

## setMethod("DatPkgFactory", "missing", function(chip) {
##   new("AffyDatPkg", name="UNKNOWN")
## })

.strMatch <- function(pat, s){length(grep(pat, s)) > 0}

setMethod("DatPkgFactory", "character", function(chip){
  if(.strMatch(".db$",chip)) chip<- sub(".db","",chip)
  pkg = paste(chip,".db",sep="")    
  ##Use standardized schema names to decide
  if(require(pkg, character.only = TRUE)){
    conn <- do.call(paste(chip, "_dbconn", sep=""), list())
    schema <- dbmeta(conn, "DBSCHEMA")    
    if(schema == "YEAST_DB" || schema == "YEASTCHIP_DB")
      pkg <- new("YeastDatPkg", name=chip)
    else if( schema == "ARABIDOPSIS_DB" || schema == "ARABIDOPSISCHIP_DB" )
      pkg <- new("ArabidopsisDatPkg", name=chip)
    else if( .strMatch("CHIP_DB$", schema)){
      pkg <- new("AffyDatPkg", name=chip)}
    else { ##Otherwise its an ordinary org package
      pkg <- new("Org.XX.egDatPkg", name=chip)
    }
    return(pkg)
  }else stop(paste("Required annotation package", chip, "is not available.",sep=" "))
})




####################################################################
## Classes and constructors to support use of GSEABase objects inside of GOstats:

setClass("GeneSetCollectionAnnotation", contains="character")

.GeneSetCollectionAnnotation <- function(annotation)
    new("GeneSetCollectionAnnotation", annotation)


GeneSetCollectionDatPkg <- function(GeneSetCollection) 
{
  GSCTypeWarning = paste("In order for the analysis to work properly, the ",
    "GeneSetCollection object must be based upon a GO2ALL mapping.  Using the ",
    "GeneSetCollection constructor that starts with a GOAllFrame will help to ",
    "guarantee this.",sep="")
  if(class(geneIdType(GeneSetCollection@.Data[[1]]))!="GOAllFrameIdentifier")
    {stop(paste(strwrap(GSCTypeWarning, exdent=2),collapse="\n"))}
    new("GeneSetCollectionDatPkg",
        GeneSetCollection=GeneSetCollection)
}


## Constructor function for parameter object needed by GOstats
GSEAGOHyperGParams <- function(name, gsc, geneIds, universeGeneIds,
                               ontology, pvalueCutoff, conditional,
                               testDirection, ...) {
    sizeWarning = paste("There is no data in your GeneSetCollection object ",
      "to do any analysis with.", sep="")
    if(length(gsc)==0){stop(paste(strwrap(sizeWarning, exdent=2),collapse="\n"))}
    new("GOHyperGParams",
        geneIds=geneIds,
        universeGeneIds=universeGeneIds,
        ontology=ontology,
        annotation=.GeneSetCollectionAnnotation(name),
        datPkg=GeneSetCollectionDatPkg(gsc),
        pvalueCutoff=pvalueCutoff,
        conditional=conditional,
        testDirection=testDirection,
        ...)
}




####################################################################
## configureDatPkg methods

setMethod("configureDatPkg", "missing", function(annotation, object) {
  object@datPkg <- DatPkgFactory(annotation)})

setMethod("configureDatPkg", "character", function(annotation, object) {
  object@datPkg <- DatPkgFactory(annotation)})

setMethod("configureDatPkg", "GeneSetCollectionAnnotation", function(annotation, object){
  if(class(object@datPkg)!= "GeneSetCollectionDatPkg"){stop("datPkg must be a GeneSetCollectionDatPkg")}
})


