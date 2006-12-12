setMethod("initialize", "HyperGParams",
          function(.Object, ...) {
              .Object <- callNextMethod(.Object, ...)
              .Object@datPkg <- DatPkgFactory(.Object@annotation)
              makeValidParams(.Object)
          })
