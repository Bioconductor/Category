setMethod("categoryName", signature(p="HyperGParams"),
          function(p) {
              p@categoryName
          })

setMethod("categoryName", signature(p="GOHyperGParams"),
          function(p) {
              c(p@categoryName, p@ontology)
          })

