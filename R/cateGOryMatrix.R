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

  if(sparse) {
    
    ## wh 14.1.2006: I considered using the matrix.csr class from the SparseM package
    ## here, but it doesn't provide row- and column names, which would make the return
    ## of 'gocats' (GO category names) and 'genes' (gene names) difficult.

    ## Really the graph is a bipartite graph, but I use graphNEL for now.
    ## The fastest way to construct it seems to be via the from-to matrix ft

    ## This code is elegant but slow:
    ##   ft = do.call("rbind", args=mapply(cbind, x, categAnc))
    ##
    ## This is a bit explicit but seems much faster:
    rg = c(0, cumsum(listLen(categAnc)))
    ft = matrix("", nrow=rg[length(rg)], ncol=2)
    for(j in 1:length(x)) {
      ft[ (rg[j]+1):rg[j+1], 1] = x[j]
      ft[ (rg[j]+1):rg[j+1], 2] = categAnc[[j]]
    }
    
    ## remove duplicated edges
    ft = ft[!duplicated(paste(ft[,1], ft[,2])), ]
    res = ftM2graphNEL(ft, edgemode="undirected")
    
  } else {
    res = matrix(as.integer(0), nrow=length(gocats), ncol=length(genes))
    rownames(res) = gocats
    colnames(res) = genes
    for(j in seq(along=x))
      res[categAnc[[j]], x[j]] = as.integer(1)
  }
  
  return(res)
}
