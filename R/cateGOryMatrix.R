## x: a vector of GO categories
augmentByAncestors = function(x) {

   s1 = x %in% ls(GOMFANCESTOR)
   s2 = x %in% ls(GOBPANCESTOR)
   s3 = x %in% ls(GOCCANCESTOR)
   s4 = is.na(x)
   
   notFound = which(!(s1|s2|s3|s4))
   if(length(notFound)>0) {
     got = if(length(notFound)>7) {
       paste(paste(notFound[1:7], collapse=", "), ", ... ")
     } else {
       paste(notFound, collapse=", ")
     }
     stop(if(length(notFound)==1) {
       sprintf("GO term %s is mentioned in 'x' but was not found in the GO package.\n", got)
     } else {
       sprintf("GO terms %s were mentioned in 'x' but were not found in the GO package.\n", got)
     })
   }
   
   stopifnot(all(s1+s2+s3+s4==1))

   res = vector(length=length(x), mode="list")
   res[s1] = mget(x[s1], GOMFANCESTOR)
   res[s2] = mget(x[s2], GOBPANCESTOR)
   res[s3] = mget(x[s3], GOCCANCESTOR)
   res[s4] = NULL
   ## the genes in x itself are not their own ancestors,
   ## so we need to add them explicitely
   res     = mapply(c, res, x) 
   names(res) = x
   return(res)
}

augmentByAncestorsSmart = function(x) {
   ## need to do the expensive search only once for each unique ID
   ux = unique(x)
   mt = match(x, ux)
   augmentByAncestors(ux)[mt]
}


## x: vector of gene identifiers (need not be unique)
## GOcat: vector of same length as x, with GO categories
cateGOry  = function(x, categ, sparse=TRUE) {

  if(!is.character(x))
    stop("'x' must be a character vector")
  if(!is.character(categ))
    stop("'categ' must be a character vector")
  if(length(x)!=length(categ))
    stop("length of 'x' and 'categ' must be the same")
  
  categAnc = augmentByAncestorsSmart(categ)

  gocats = sort(unique(unlist(categAnc)))
  genes  = sort(unique(unlist(x)))

  res = do.call(if(sparse) "Matrix" else "matrix",
    list(FALSE, nrow=length(gocats), ncol=length(genes)))
  
  rownames(res) = gocats
  colnames(res) = genes
  for(j in seq(along=x))
    res[categAnc[[j]], x[j]] = TRUE
  
  return(res)
}
