setMethod("initialize", "HyperGParams",
          function(.Object, ...) {
              .Object <- callNextMethod(.Object, ...)
              configureDatPkg(annotation(.Object), .Object)
              makeValidParams(.Object)
          })

