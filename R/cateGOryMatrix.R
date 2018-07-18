## x: a character vector of GO category IDs
augmentByAncestors = function(x) {
    if (!requireNamespace("GO.db"))
        stop("Use 'BiocManager::install(\"GO.db\")' to install the GO.db package")
   ## NOTE: this function implicitly does require("GO.db")
   s1 = x %in% ls(getAnnMap("MFANCESTOR", "GO", load=TRUE))
   s2 = x %in% ls(getAnnMap("BPANCESTOR", "GO"))
   s3 = x %in% ls(getAnnMap("CCANCESTOR", "GO"))
   s4 = is.na(x)

   ## throw warning if some not found
   notFound = x[!(s1|s2|s3|s4)]
   if(length(notFound)>0) {
     notFound[notFound==""] = "(zero length string)"
     got = if(length(notFound)>7) {
       paste(paste(notFound[1:7], collapse=", "), ", ... ")
     } else {
       paste(notFound, collapse=", ")
     }
     warning(if(length(notFound)==1) {
       sprintf("cateGOry: GO term %s is mentioned in 'x' but was not found in the GO package.\n", got)
     } else {
       sprintf("cateGOry: GO terms %s were mentioned in 'x' but were not found in the GO package.\n", got)
     })
   }

   res = vector(length=length(x), mode="list")
   res[s1] = mget(x[s1], GO.db::GOMFANCESTOR)
   res[s2] = mget(x[s2], GO.db::GOBPANCESTOR)
   res[s3] = mget(x[s3], GO.db::GOCCANCESTOR)
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
cateGOry  = function(x, categ, sparse=FALSE) {

  if(!is.character(x))
    stop("'x' must be a character vector")
  if(!is.character(categ))
    stop("'categ' must be a character vector")
  if(length(x)!=length(categ))
    stop("length of 'x' and 'categ' must be the same")

  categAnc = augmentByAncestorsSmart(categ)

  gocats = sort(unique(unlist(categAnc)))
  genes  = sort(unique(unlist(x)))

  res = do.call(if(sparse) Matrix::Matrix else base::matrix,
    list(FALSE, nrow=length(gocats), ncol=length(genes)))

  ## res = matrix(FALSE, nrow=length(gocats), ncol=length(genes))

  rownames(res) = gocats
  colnames(res) = genes
  for(j in seq(along=x))
    res[categAnc[[j]], x[j]] = TRUE

  return(res)
}
