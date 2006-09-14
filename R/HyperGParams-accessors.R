setMethod("geneIds", "HyperGParams", function(r) r@geneIds)
setReplaceMethod("geneIds", "HyperGParams", function(r, value) {
    r@geneIds <- value
    r
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
setReplaceMethod("annotation", "HyperGParams", function(object, value) {
    object@annotation <- value
    object
})

setMethod("pvalueCutoff", "HyperGParams", function(r) r@pvalueCutoff)
setReplaceMethod("pvalueCutoff", "HyperGParams", function(r, value) {
    r@pvalueCutoff <- value
    r
})


setMethod("isConditional", "HyperGParams", function(r) FALSE)


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

