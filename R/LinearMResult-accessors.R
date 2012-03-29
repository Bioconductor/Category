### Accessor methods for LinearMResultBase class

## generic "annotation" defined in Biobase

setMethod("annotation", signature(object="LinearMResultBase"),
          function(object) object@annotation)

setMethod("universeCounts", signature(r="LinearMResultBase"),
          function(r) {
              univ <- geneIdUniverse(r)
              ans <- listLen(univ)
              names(ans) <- names(univ)
              ans
          })

setMethod("universeMappedCount", signature(r="LinearMResultBase"),
           function(r) length(unique(unlist(geneIdUniverse(r)))))

setMethod("geneMappedCount", signature(r="LinearMResultBase"),
           function(r) length(geneIds(r)))

setMethod("geneIds", signature(object="LinearMResultBase"),
          function(object, ...) object@geneIds)

setMethod("geneIdsByCategory", signature(r="LinearMResultBase"),
          function(r, catids=NULL) {
              ans <- geneIdUniverse(r)
              if (!missing(catids) && !is.null(catids))
                ans <- ans[catids]
              lapply(ans, intersect, geneIds(r))
          })

setMethod("sigCategories", signature(r="LinearMResultBase"),
          function(r, p) {
              if (missing(p))
                  p <- pvalueCutoff(r)
              pv <- pvalues(r)
              names(pv[pv < p])
          })

setMethod("testName", signature(r="LinearMResultBase"),
          function(r) r@testName)

setMethod("pvalueCutoff", signature(r="LinearMResultBase"),
          function(r) r@pvalueCutoff)

setMethod("testDirection", signature(r="LinearMResultBase"),
          function(r) r@testDirection)

setMethod("description",
          signature(object="LinearMResultBase"),
          function(object) {
              cond <- "Conditional"
              if (!conditional(object))
                cond <- ""
              desc <- paste("Gene to %s ", cond, " test for %s-regulation", sep = "")
              desc <- sprintf(desc, paste(testName(object), collapse=" "),
                              testDirection(object))
              desc
          })

setMethod("conditional", "LinearMResultBase",
          function(r) {
              if (!("conditional" %in% slotNames(r)))
                  FALSE
              else
                  r@conditional
          })

setMethod("isConditional", "LinearMResultBase",
          function(r) conditional(r))


setMethod("condGeneIdUniverse", signature(r="LinearMResultBase"),
          function(r) geneIdUniverse(r, cond=TRUE))

setMethod("geneIdUniverse", signature(r="LinearMResultBase"),
          function(r, cond=TRUE) {
            im <- addHierarchyToIncidenceMatrix(incidence(r@gsc), r@graph)
            split(colnames(im)[col(im)[as.logical(im)]],
                  factor(rownames(im)[row(im)[as.logical(im)]], nodes(r@graph)))
          })

### Accessor methods for LinearMResult class

setMethod("pvalues", signature(r="LinearMResult"),
          function(r) r@pvalues)

setMethod("effectSize", signature(r="LinearMResult"),
          function(r) r@effectSize)
