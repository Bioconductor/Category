.geneSetParamListFlip <- function(p){
    ## yes, it is already a list, but its backwards of what we
    ## want here.
    coll <- p@geneSetCollection
    genes <- geneIds(coll)
    genesLengths <- lapply(genes, length)
    IDs <- names(coll)
    IDReps <- rep(IDs, genesLengths)
    collFrame <- cbind(IDReps, unlist(genes))
    collList <- split(as.character(collFrame[,1]),
                      as.character(collFrame[,2]))
    return(collList)
}

getDbAnnMap <- function(col, db){
    if(col %in% columns(db))
        return(AnnotationDbi:::makeFlatBimapUsingSelect(db, col))
    else
        stop(paste("The database", db, "you are using has no", col, "data."), call. = FALSE)
}

annMapChooser <- function(col, datpkg){
    if(datpkg@installed)
        getAnnMap(col, datpkg@name)
    else
        getDbAnnMap(col, datpkg@db)
}

setMethod("ID2GO", "DatPkg",
          function(p) annMapChooser("GO", p))

setMethod("ID2GO", "GeneSetCollectionDatPkg",
          function(p) list2env(.geneSetParamListFlip(p)))

setMethod("ID2KEGG", "DatPkg",
          function(p) annMapChooser("PATH", p))

setMethod("ID2KEGG", "GeneSetCollectionDatPkg",
          function(p) list2env(.geneSetParamListFlip(p)))




setMethod("ID2EntrezID", "AffyDatPkg",
          function(p) annMapChooser("ENTREZID", p))


##FIXME: this is seriously slow - try list2env to speed up a bit
.createIdentityMap <- function(keys) {
    keys = as.list(keys)
    names(keys) = keys 
    list2env(keys)
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
              else if(p@installed)
                  return(.createIdentityMap(allValidKeys(p@name)))
              else
                  .createIdentityMap(allValidKeys(p@db))
})

setMethod("ID2EntrezID", "ArabidopsisDatPkg",
          function(p) {
              bname = p@name 
              if( exists( paste(bname, "ACCNUM", sep="")) ) 
	        return(getAnnMap("ACCNUM", p@name))
              else if(p@installed)
                  return(.createIdentityMap(allValidKeys(p@name)))
              else
                  .createIdentityMap(allValidKeys(p@db))
})

setMethod("ID2EntrezID", "Org.XX.egDatPkg",
          function(p) {
    if(p@installed)
        return(.createIdentityMap(allValidKeys(p@name)))
    else
        .createIdentityMap(allValidKeys(p@db))
})


setMethod("ID2EntrezID", "GeneSetCollectionDatPkg", function(p) {
    ## This method does not need to really do anything "real" since
    ## they are going to get out the ID type that they put in: no
    ## matter what.
    coll <- p@geneSetCollection
    genes <- unique(unlist(geneIds(coll)))
    collList <- split(genes,genes)
    res <- list2env(collList)
    res
})

setMethod("GO2AllProbes", "DatPkg",
          function(p, ontology=c("BP", "CC", "MF")) {
    ontIds <- aqListGOIDs(ontology)
    if(p@installed)
        go2all <- getAnnMap("GO2ALLPROBES", p@name)
    else
        go2all <- getDbAnnMap("GOALL", p@db)
    ontIds <- intersect(ontIds, ls(go2all))
    go2allOnt <- mget(ontIds, go2all, ifnotfound=NA)
    go2allOnt <- removeLengthZeroAndMissing(go2allOnt)
    list2env(go2allOnt)
})


setMethod("GO2AllProbes", "YeastDatPkg",
          function(p, ontology=c("BP", "CC", "MF")) {
    conn <- dbconn(p@db)
    schema <- dbmeta(conn, "DBSCHEMA")
    env = environment()
    if(schema == "YEASTCHIP_DB"){
        env = callNextMethod()
        return(env)
    }
    ontIds <- aqListGOIDs(ontology)
    if(p@installed)
        go2all <- getAnnMap("GO2ALLORFS", p@name)
    else
        go2all <- getDbAnnMap("GOALL", p@db)
    ontIds <- intersect(ontIds, ls(go2all))
    go2allOnt <- mget(ontIds, go2all, ifnotfound=NA)
    go2allOnt <- removeLengthZeroAndMissing(go2allOnt)
    env = list2env(go2allOnt)
    return(env)
})





setMethod("GO2AllProbes", "Org.XX.egDatPkg",
          function(p, ontology=c("BP", "CC", "MF")) {

              #db <- get("db_conn", paste("package:", p@name, sep=""))
              db <- dbconn(p@db)
              sqlQ <- "SELECT DISTINCT gene_id, go_id
              FROM genes INNER JOIN go_%s USING (_id)"
              sqlQ <- sprintf(sqlQ, tolower(ontology))
              go2egTab <- dbGetQuery(db, sqlQ)
              go2eg <- list2env(split(go2egTab[["gene_id"]], go2egTab[["go_id"]]))

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
            coll <- p@geneSetCollection
            ## Lets put the GeneSetCollection into a format that is
            ## easier to filter (from the left OR right)
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
            result <- split(as.character(collFrame[,2]),
                            as.character(collFrame[,1]))
            
            if(length(result)==0)
              stop("no annotations for selected genes")
            list2env(result)
          })




setMethod("KEGG2AllProbes", "DatPkg",
          function(p) revmap(annMapChooser("PATH", p)))

setMethod("KEGG2AllProbes", "GeneSetCollectionDatPkg",  
          function(p) {
            coll <- p@geneSetCollection
            ## Lets put the GeneSetCollection into a format that is
            ## easier to filter (from the left OR right)
            genes = geneIds(coll)
            genesLengths = lapply(genes, length)
            KEGGIDs = names(coll)
            KEGGIDReps = rep(KEGGIDs, genesLengths)
            collFrame = cbind(KEGGIDReps, unlist(genes))
            
            ##Then put things back into a list format
            result <- split(as.character(collFrame[,2]),
                            as.character(collFrame[,1]))
            
            if(length(result)==0)
              stop("no annotations for selected genes")
            list2env(result)
          })



setMethod("isDBDatPkg","DatPkg",
          function(d) length(d@db) > 0)

















####################################################################
## DatPkgFactory methods (default has to be back in AllClasses) :(

## setMethod("DatPkgFactory", "missing", function(chip) {
##   new("AffyDatPkg", name="UNKNOWN")
## })

.strMatch <- function(pat, s){length(grep(pat, s)) > 0}

setMethod("DatPkgFactory", "character", function(chip) {
    if (.strMatch(".db$",chip))
        chip<- sub(".db","",chip)
    pkg <- paste(chip,".db",sep="")    
    if(!require(pkg, character.only = TRUE))
        stop("annotation package '", pkg, "' not available")
    
    ## Use standardized schema names to decide
    db <- get(pkg)
    conn <- dbconn(db)
    schema <- dbmeta(conn, "DBSCHEMA")
    if (schema == "YEAST_DB" || schema == "YEASTCHIP_DB")
        new("YeastDatPkg", name=chip, db=db, installed=TRUE)
    else if( schema == "ARABIDOPSIS_DB" || schema == "ARABIDOPSISCHIP_DB" )
        new("ArabidopsisDatPkg", name=chip, db=db, installed=TRUE)
    else if( .strMatch("CHIP_DB$", schema))
        new("AffyDatPkg", name=chip, db=db, installed=TRUE)
    else ## Otherwise its an ordinary org package
        new("Org.XX.egDatPkg", name=chip, db=db, installed=TRUE)
})

setMethod("DatPkgFactory", "OrgDb", function(chip) {
    ## don't act like installed package is from AnnotationHub
    pkName <- get("packageName", chip@.xData)
    if(length(pkName) > 0)
        return(DatPkgFactory(pkName))
    conn <- dbconn(chip)
    schema <- dbmeta(conn, "DBSCHEMA")
    orgn <- dbmeta(conn, "ORGANISM")
    chiptype <- dbmeta(conn, "Db type")
    name <- paste(chiptype, "for", orgn)
    if (schema == "YEAST_DB")
        new("YeastDatPkg", name=name, db=chip, installed=FALSE)
    else if( schema == "ARABIDOPSIS_DB")
        new("ArabidopsisDatPkg", name=name, db=chip, installed=FALSE)
    else ## Otherwise its an ordinary org package
        new("Org.XX.egDatPkg", name=name, db=chip, installed=FALSE)
})

setMethod("DatPkgFactory", "ChipDb", function(chip) {
    ## don't act like installed package is from AnnotationHub
    pkName <- get("packageName", chip@.xData)
    if(length(pkName) > 0)
        return(DatPkgFactory(pkName))
    conn <- dbconn(chip)
    schema <- dbmeta(conn, "DBSCHEMA")
    orgn <- dbmeta(conn, "ORGANISM")
    chiptype <- dbmeta(conn, "Db type")
    name <- paste(chiptype, "for", orgn)
    if (schema == "YEASTCHIP_DB")
        new("YeastDatPkg", name=name, db=chip, installed=FALSE)
    else if(schema == "ARABIDOPSISCHIP_DB" )
        new("ArabidopsisDatPkg", name=name, db=chip, installed=FALSE)
    else
        new("AffyDatPkg", name=name, db=chip, installed=FALSE)
})



## this is for OBO
OBOCollectionDatPkg <- function(oboCollection, geneSetCollection) {
  new("OBOCollectionDatPkg", oboCollection=oboCollection,
      oboGraph=as(oboCollection, "graphNEL"),
      geneSetCollection=geneSetCollection,
      installed = FALSE)
}
                                       
OBOHyperGParams <- function(...)
    new("OBOHyperGParams", ...)

####################################################################
## Classes and constructors to support use of GSEABase objects inside
## of GOstats:

setClass("GeneSetCollectionAnnotation", contains="character")

.GeneSetCollectionAnnotation <- function(annotation)
    new("GeneSetCollectionAnnotation", annotation)


GeneSetCollectionDatPkg <- function(geneSetCollection) 
{
    new("GeneSetCollectionDatPkg",
        geneSetCollection=geneSetCollection, installed = FALSE)
}




## Constructor function for parameter object needed by GOstats
GSEAGOHyperGParams <-
    function(name, geneSetCollection, geneIds, universeGeneIds,
             ontology, pvalueCutoff, conditional, testDirection, ...)
{
    if(!extends(class(geneIdType(geneSetCollection[[1]])),
                "GOAllFrameIdentifier"))
    {
        GSCTypeWarning <-
            paste("'geneSetCollection' elements must use GO2ALL",
                  "mappings; use GeneSetCollection constructors",
                  "that start with 'GOAllFrame'")
        stop(paste(strwrap(GSCTypeWarning, exdent=2),collapse="\n"))
    }
    if(length(geneSetCollection)==0)
        stop("geneSetCollection has length 0")
    new("GOHyperGParams",
        geneIds=geneIds,
        universeGeneIds=universeGeneIds,
        ontology=ontology,
        annotation=.GeneSetCollectionAnnotation(name),
        datPkg=GeneSetCollectionDatPkg(geneSetCollection),
        pvalueCutoff=pvalueCutoff,
        conditional=conditional,
        testDirection=testDirection,
        ...)
}

## Constructor function for parameter object needed by GOstats
GSEAKEGGHyperGParams <-
    function(name, geneSetCollection, geneIds, universeGeneIds,
             pvalueCutoff, testDirection, ...)
{

    if(length(geneSetCollection)==0)
        stop("geneSetCollection has length 0")
    new("KEGGHyperGParams",
        geneIds=geneIds,
        universeGeneIds=universeGeneIds,
        annotation=.GeneSetCollectionAnnotation(name),
        datPkg=GeneSetCollectionDatPkg(geneSetCollection),
        pvalueCutoff=pvalueCutoff,
        testDirection=testDirection,
        ...)
}


####################################################################
## configureDatPkg methods

setMethod("configureDatPkg", "character",
          function(annotation, ...) DatPkgFactory(annotation))

setMethod("configureDatPkg", "GeneSetCollectionAnnotation",
          function(annotation, object, ...) object@datPkg)

setMethod("configureDatPkg", "OrgDb",
          function(annotation, ...) DatPkgFactory(annotation))



####################################################################
## organism method

setMethod("organism", "HyperGParams",
          function(object) organism(object@datPkg) )


setMethod("organism", "GeneSetCollectionDatPkg",
          function(object) organism(object@geneSetCollection[[1]]) )

setMethod("organism", "DatPkg",
          function(object)
    if(object@installed) organism(object@name) else dbmeta(dbconn(object@db), "ORGANISM"))


setMethod("organism", "HyperGResult",
          function(object) object@organism )
