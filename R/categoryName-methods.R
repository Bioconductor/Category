setMethod("categoryName", signature(p="GeneCategoryHyperGeoTestParams"),
          function(p) {
              p@categoryName
          })

setMethod("categoryName", signature(p="GeneGoHyperGeoTestParams"),
          function(p) {
              c(p@categoryName, p@ontology)
          })

