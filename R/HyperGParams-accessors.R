## autogenerate accessors

if (FALSE) {
theSlots <- slotNames("HyperGParams")
for (s in theSlots) {
    if (!isGeneric(s))
      setGeneric(s, function(r) standardGeneric(s))
    if (!isGeneric(paste(s, "<-", sep="")))
      setGeneric(paste(s, "<-", sep=""),
                 function(r, value) standardGeneric(paste(s, "<-", sep="")))
    
    setMethod(s, signature(r="HyperGParams"),
              function(r) slot(r, s))
    
    setReplaceMethod(s,
                     signature(r="HyperGParams"),
                     function(r, value) {
                         slot(r, s) <- value
                         r
                     })
}
              
}
