setMethod("initialize", "HyperGParams",
          function(.Object, ...)
{
    .Object <- callNextMethod()
    .Object@datPkg <- 
        configureDatPkg(annotation(.Object), .Object)
    makeValidParams(.Object)
})

