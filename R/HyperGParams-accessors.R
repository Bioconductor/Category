.makeValidParams <- function(object) {
    ## check HyperGParams instance for validity.
    ## If we can fix it, we do (and issue a warning)
    ## Return a more valid instance or error

    ##Check if annotation has been written "long form"
    ##If it is, then shorten the name appropriately.
    ann <- annotation(object)
    if(length(grep(".db$", ann)) > 0){
        ann<- sub("\\.db$", "", ann)
        annotation(object) <- ann
    }
    
    sel <- geneIds(object)
    if (is.list(sel)) {
        warning("converting geneIds from list to atomic vector via unlist")
        sel <- unlist(sel)
    }
    if (any(duplicated(sel))) {
        warning("removing duplicate IDs in geneIds")
        sel <- unique(sel)
    }
    geneIds(object) <- sel
    univ <- universeGeneIds(object)
    if (length(univ)) {
        if (is.list(univ)) {
            warning("converting univ from list to atomic vector via unlist")
            univ <- unlist(univ)
        }
        if (typeof(sel) != typeof(univ))
          stop(paste("geneIds and universeGeneIds must have the same mode\n",
                       "geneIds:", typeof(sel), "\n",
                       "universeGeneIds:", typeof(univ)), .Call=FALSE)
        if (any(duplicated(univ))) {
            warning("removing duplicate IDs in universeGeneIds")
            univ <- unique(univ)
        }
        universeGeneIds(object) <- univ
        if (!all(sel %in% univ)) {
            warning("removing geneIds not in universeGeneIds")
            sel <- intersect(sel, univ)
            if (!length(sel))
              stop("no geneIds in universeGeneIds", .Call=FALSE)
            geneIds(object) <- sel
        }
    }
    pv <- pvalueCutoff(object)
    if (pv > 1 || pv < 0)
      stop("invalid pvalueCutoff, must be between 0 and 1", .Call=FALSE)
    if (length(annotation(object)) != 1)
      stop("annotation must be a length 1 character vector", .Call=FALSE)
    object
}
setMethod("makeValidParams", "HyperGParams", .makeValidParams)

setMethod("geneIds", "HyperGParams", function(object, ...) object@geneIds)
setReplaceMethod("geneIds", "HyperGParams", function(object, value) {
    object@geneIds <- value
    object
})

setMethod("categorySubsetIds", "HyperGParams", function(r) r@categorySubsetIds)
setReplaceMethod("categorySubsetIds", "HyperGParams", function(r, value) {
    r@categorySubsetIds <- value
    r
})

setMethod("testDirection", "HyperGParams", function(r) r@testDirection)
setReplaceMethod("testDirection", "HyperGParams", function(r, value) {
    r@testDirection <- value
    r
})

setMethod("universeGeneIds", "HyperGParams", function(r) r@universeGeneIds)
setReplaceMethod("universeGeneIds", "HyperGParams", function(r, value) {
    r@universeGeneIds <- value
    r
})

setMethod("categoryName", "HyperGParams", function(r) r@categoryName)
setReplaceMethod("categoryName", "HyperGParams", function(r, value) {
    r@categoryName <- value
    r
})

setMethod("annotation", "HyperGParams", function(object) object@annotation)
setReplaceMethod("annotation", c("HyperGParams", "character"),
                 function(object, value) {
                   object@annotation <- value
                   object@datPkg <- DatPkgFactory(value)
                   object
                 })

setMethod("pvalueCutoff", "HyperGParams", function(r) r@pvalueCutoff)
setReplaceMethod("pvalueCutoff", "HyperGParams", function(r, value) {
    r@pvalueCutoff <- value
    r
})


setMethod("conditional", "HyperGParams", function(r) FALSE)

setMethod("conditional", "ChrMapHyperGParams", function(r) r@conditional)

setReplaceMethod("conditional", c("ChrMapHyperGParams", "logical"),
                 function(r, value) {
                     if (is.na(value))
                       stop("value must be TRUE or FALSE")
                     r@conditional <- value
                     r
                 })

setMethod("conditional", "GOHyperGParams", function(r) r@conditional)

setReplaceMethod("conditional", c("GOHyperGParams", "logical"),
                 function(r, value) {
                     if (is.na(value))
                       stop("value must be TRUE or FALSE")
                     r@conditional <- value
                     r
                 })

setMethod("isConditional", "GOHyperGParams", function(r) conditional(r))

setMethod("ontology", "HyperGParams", function(object) NA)

setMethod("ontology", "GOHyperGParams", function(object) object@ontology)

setReplaceMethod("ontology", c("GOHyperGParams", "character"),
                 function(r, value) {
                     if (is.na(value) || length(value) != 1)
                       stop("value must be a length one character vector")
                     r@ontology <- value
                     r
                 })

##FIXME, this shouldn't be as hard as it is :-( :-(
## autogenerate accessors
## theSlots <- slotNames("HyperGParams")

## getter <- function(r) slot(r, s)
## getterArgs <- alist(r=)

## setter <- function(r, v) {
##     slot(r, s) <- v
##     r
## }
## setterArgs <- alist(r=, v=)

## for (s in theSlots) {
##     cat("Creating get method for ", s, "\n")
##     newGetArgs <- eval(substitute(c(r),
##                                   list(r=getGeneric(s)@signature[1])))
##     names(getterArgs) <- newGetArgs
##     formals(getter) <- getterArgs
##     setMethod(s, signature("HyperGParams"), getter)

##     sSET <- paste(s, "<-", sep="")
##     cat("Creating set method for ", s, "\n")
##     newSetArgs <- eval(substitute(c(r, v),
##                                   list(r=getGeneric(sSET)@signature[1],
##                                        v=getGeneric(sSET)@signature[2])))
##     names(setterArgs) <- newSetArgs
##     formals(setter) <- setterArgs
    
##     setReplaceMethod(s, signature("HyperGParams"), setter)
## }

