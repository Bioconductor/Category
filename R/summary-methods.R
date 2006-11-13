getWantedResults <- function(result, pvalue, categorySize=NULL) {
    ## Returns a logical vector with TRUE indicating selected
    ## results from those tested in the specified result instance.
    pvals <- pvalues(result)
    wanted <- pvals < pvalue
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


htmlReportFromDf <- function(r, caption, file="", append=FALSE)
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
    dig <- rep(2, ncol(r)+1)  ## need +1 for xtable row name
    dig[5:6] <- 0
    xt <- xtable(r, caption=caption,
                 digits=dig)
    print(xt, type="html", file=file, append=append,
          caption.placement="top")
    return(invisible(TRUE))
}


setMethod("htmlReport", signature(r="HyperGResultBase"),
          function(r, file="", append=FALSE, label="", ...)
          {
              htmlReportFromDf(r=summary(r, ...),
                               caption=paste(label, description(r)),
                               file=file, append=append)
          })

