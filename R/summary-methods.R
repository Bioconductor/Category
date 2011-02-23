getWantedResults <-
    function(result, pvalue, categorySize=NULL)
{
    ## Returns a logical vector with TRUE indicating selected
    ## results from those tested in the specified result instance.
    pvals <- pvalues(result)
    wanted <- is.finite(pvals) & pvals < pvalue
    if (!is.null(categorySize)) {
        ucounts <- universeCounts(result)
        hasMinSize <- ucounts >= categorySize
        wanted <- wanted & hasMinSize
    }
    wanted
}

setMethod("summary", signature(object="HyperGResultBase"),
          function(object, pvalue=pvalueCutoff(object), categorySize=NULL)
          {
              ## Filter based on p-value and category size
              wanted <- getWantedResults(object, pvalue, categorySize)
              pvals <- pvalues(object)
              ucounts <- universeCounts(object)
              if (!any(wanted)) {
                  warning("No results met the specified criteria.  ",
                          "Returning 0-row data.frame", call.=FALSE)
                  catIds <- character(0)
                  pvals <- odds <- ecounts <- numeric(0)
                  counts <- ucounts <- integer(0)
              } else {
                  pvals <- pvals[wanted]
                  ucounts <- ucounts[wanted]
                  catIds <- names(pvals)
                  odds <- oddsRatios(object)[wanted]
                  ecounts <- expectedCounts(object)[wanted]
                  counts <- geneCounts(object)[wanted]
              }
              df <- data.frame(ID=catIds, Pvalue=pvals, OddsRatio=odds,
                               ExpCount=ecounts, Count=counts, Size=ucounts,
                               stringsAsFactors=FALSE, row.names=NULL)
              names(df)[1] <- paste(paste(testName(object), collapse=""),
                                    "ID", sep="")
              df
          })

setMethod("summary", signature(object="KEGGHyperGResult"),
          function(object, pvalue=pvalueCutoff(object),
                   categorySize=NULL, htmlLinks=FALSE){
              KEGG_URL <- "http://www.genome.jp/dbget-bin/www_bget?path:%s%s"
              ## annOrg <- get(paste(annotation(object), "ORGANISM", sep=""))
              annOrg <- organism(object)
              orgSpecifier <- switch(annOrg,
                                     "Homo sapiens"="hsa",
                                     "Mus musculus"="mmu",
                                     "Rattus norvegicus"="rnu",
                                     ## will need others in future
                                     "hsa")
              df <- callNextMethod(object=object, pvalue=pvalue,
                                   categorySize=categorySize)
              if(nrow(df) == 0){
                  df$Term <- character(0)
                  return(df)
              }
              keggIds <- df[[1]]
              ## implicit require("KEGG.db")
              keggEnv <- getAnnMap("PATHID2NAME", "KEGG", load=TRUE)
              keggTerms <- unlist(mget(keggIds, keggEnv, ifnotfound=NA))
              if(htmlLinks){
                  keggIdUrls <- sapply(keggIds,
                                       function(x)
                                       sprintf(KEGG_URL, orgSpecifier, x))
                  keggTerms <- paste('<a href="', keggIdUrls, '">', keggTerms,
                                     '</a>', sep="")
              }
              df$Term <- keggTerms
              df
          })

setMethod("summary", signature(object="PFAMHyperGResult"),
          function(object,pvalue=pvalueCutoff(object),
                   categorySize=NULL, htmlLinks=FALSE){
              PFAM_URL <- "http://pfam.sanger.ac.uk/family?acc=%s"
              df <- callNextMethod(object=object, pvalue=pvalue,
                                   categorySize=categorySize)
              if(nrow(df) == 0){
                  df$Term <- character(0)
                  return(df)
              }
              pfamIds <- df[[1]]
              if(htmlLinks){
                  pfamIdUrls <- sapply(pfamIds,
                                       function(x) sprintf(PFAM_URL, x))
                  pfamTerms <- paste('<a href="', pfamIdUrls, '">', pfamIds,
                                     '</a>', sep="")
              }else{
                  pfamTerms <- pfamIds
              }
              df$Term <- pfamTerms
              df
          })


htmlReportFromDf <- function(r, caption, file="", append=FALSE, digits=3)
{
    have_xtable <- suppressWarnings({
        require("xtable", quietly=TRUE, warn.conflicts=FALSE)
    })
    if (!have_xtable)
      stop("htmlReport needs the xtable package and it is not",
           "available", call.=FALSE)
    if (nrow(r) == 0) {
        warning("No rows to report.  Skipping")
        return(invisible(FALSE))
    }
    ## XXX: Hard-coded column formatting here
    dig <- rep(digits, ncol(r)+1)  ## need +1 for xtable row name
    dig[5:7] <- 0
    xt <- xtable(r, caption=caption,
                 digits=dig)
    print(xt, type="html", file=file, append=append,
          caption.placement="top",
          sanitize.text.function=function(x) x,
          include.rownames=FALSE)
    return(invisible(TRUE))
}

XXX_getSummaryGeneric_XXX <- function() {
    ## FIXME: the methods packge is broken for this case
    ## so we have to find the right summary method ourselves
    places <- find("summary")
    ## take the first standardGeneric
    f <- NULL
    for (i in seq(along=places)) {
        f <- get("summary", places[i])
        if (is(f, "standardGeneric"))
          break
        else
          f <- NULL
    }
    if (is.null(f))
      stop("could not find appropriate summary method")
    f
}

setMethod("htmlReport", signature(r="HyperGResultBase"),
          function(r, file="", append=FALSE, label="", digits=3,
                   summary.args=NULL)
          {
              summary <- XXX_getSummaryGeneric_XXX()
              if (!is.null(summary.args) && !is.list(summary.args))
                stop("'summary.args' must be NULL or a list of arguments for",
                     " the summary method")
              df <- do.call(summary, c(list(r), summary.args))
              htmlReportFromDf(r=df,
                               caption=paste(label, description(r)),
                               file=file, append=append, digits=digits)
          })

setMethod("htmlReport", signature(r="KEGGHyperGResult"),
          function(r, file="", append=FALSE, label="",
                   digits=3, summary.args=list(htmlLinks=TRUE)){
              callNextMethod(r=r, file=file, append=append,
                             label=label, digits=digits,
                             summary.args=summary.args)
          })

setMethod("htmlReport", signature(r="PFAMHyperGResult"),
          function(r, file="", append=FALSE, label="",
                   digits=3, summary.args=list(htmlLinks=TRUE)){
              callNextMethod(r=r, file=file, append=append,
                             label=label, digits=digits,
                             summary.args=summary.args)
          })



setMethod("summary", signature(object="LinearMResultBase"),
          function(object,
                   pvalue = pvalueCutoff(object),
                   categorySize = NULL,
                   ...)
      {
          ##               ## FIXME: should do this in a better way
          ##               object@pvalues <- p.adjust(pvalues(object), method = adjust.pvalues)

          ## Filter based on p-value and category size
          wanted <- getWantedResults(object, pvalue, categorySize)
          pvals <- pvalues(object)
          esize <- effectSize(object)
          ucounts <- universeCounts(object)
          if (!any(wanted)) {
              warning("No results met the specified criteria.  ",
                      "Returning 0-row data.frame", call.=FALSE)
              catIds <- character(0)
              pvals <- numeric(0)
              esize <- numeric(0)
              ucounts <- integer(0)
          } else {
              pvals <- pvals[wanted]
              esize <- esize[wanted]
              ucounts <- ucounts[wanted]
              catIds <- names(pvals)
          }
          df <- data.frame(ID=catIds, Pvalue=pvals,
                           Effect=esize,
                           Size=ucounts,
                           stringsAsFactors=FALSE, row.names=NULL)
          names(df)[1] <- paste(paste(testName(object), collapse=""),
                                "ID", sep="")
          df
      })

