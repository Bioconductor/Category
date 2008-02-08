setMethod("categoryName", signature(r="HyperGParams"),
          function(r) {
              r@categoryName
          })

setMethod("categoryName", signature(r="GOHyperGParams"),
          function(r) {
              c(r@categoryName, r@ontology)
          })



setMethod("categoryName", signature(r="LinearMParams"),
          function(r) {
              r@categoryName
          })
