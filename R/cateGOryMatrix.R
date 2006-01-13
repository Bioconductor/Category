## x: a vector of GO categories
augmentByAncestors = function(x) {

   s1 = x %in% ls(GOMFANCESTOR)
   s2 = x %in% ls(GOBPANCESTOR)
   s3 = x %in% ls(GOCCANCESTOR)
   s4 = is.na(x)
   
   notFound = which(!(s1|s2|s3|s4))
   if(length(notFound)>0)
     stop(sprintf("GO term%s %s not found.\n", ifelse(length(notFound)>1, "s", ""),
                  paste(x[notFound], collapse=", ")))
   
   stopifnot(all(s1+s2+s3+s4==1))

   res = vector(length=length(x), mode="list")
   res[s1] = mget(x[s1], GOMFANCESTOR)
   res[s2] = mget(x[s2], GOBPANCESTOR)
   res[s3] = mget(x[s3], GOCCANCESTOR)
   res[s4] = NA
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
cateGOryMatrix  = function(x, categ) {
  
  if(!require("GO"))
    stop("Cannot load the 'GO' package")

  categAnc = augmentByAncestorsSmart(categ)

  allGO = sort(unique(unlist(categAnc)))
  allx  = sort(unique(x))
  
  res = matrix(as.integer(0), nrow=length(allGO), ncol=length(allx))
  rownames(res) = allGO
  colnames(res) = allx
  
  for(j in seq(along=x))
    res[categAnc[[j]], x[j]] = as.integer(1)

  return(res)
}
