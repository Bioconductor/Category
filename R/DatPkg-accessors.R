setMethod("ID2GO", "DatPkg",
          function(p) getDataEnv("GO", p@name))

setMethod("GO2AllProbes", "DatPkg",
          function(p) getDataEnv("GO2ALLPROBES", p@name))

setMethod("ID2EntrezID", "AffyDatPkg",
          function(p) getDataEnv("ENTREZID", p@name))

setMethod("ID2EntrezID", "YeastDatPkg",
          function(p) {
              ## Create an identity map
              e <- new.env(parent=emptyenv(), hash=TRUE)
              for (n in ls(getDataEnv("CHR", p@name))) {
                  e[[n]] <- n
              }
              e
          })
                   
setMethod("GO2AllProbes", "OrganismMappingDatPkg",
          function(p) {
              ## Not sure, I think we have to build it up :-(
              getDataEnv("GO2ALLPROBES", p@name)
          })
              
