local_test_factory <- function(selids, PCUT=0.05, tableTest=chisq.test)
{
    force(selids)
    anyexist <- function(keys, env) {
        ans <- mget(keys, env, ifnotfound=list(NULL))
        any(sapply(ans, function(x) !is.null(x)))
    }
    
    function(start, g, prev_ans) {
        seen <- new.env(hash=TRUE, parent=emptyenv())
        for (a in prev_ans) {
            sapply(a[["nodes"]], function(i) seen[[i]] <- 1L)
        }
        ans <- lapply(start, function(sibs) {
            if (!length(sibs) || anyexist(sibs, seen)) {
                return(NULL)
            }
            sapply(sibs, function(s) seen[[s]] = 1L)
            band2gene <- lgeneIds(g, sibs)
            imat <- ann_list_to_mat(band2gene)
            ctab <- cb_buildContinTable(imat, selids, min.expected=NULL)
            if (length(ctab))
              c(list(nodes=sibs), tableTest(ctab))
            else
              list(nodes=sibs,
                   statistic=NA, parameter=NA, p.value=NA,
                   method="Test NOT PERFORMED",
                   data.name=sibs,
                   observed=NA, expected=NA, residuals=NA)
        })
        ans[!sapply(ans, is.null)]
    }
}

make_hyperg_result <- function(node, numW, numB, numDrawn, numWdrawn, OVER)
{
    n11 <- numWdrawn
    n21 <- numW - numWdrawn
    n12 <- numDrawn - numWdrawn
    n22 <- numB - n12

    obs <- matrix(c(n11, n21, n12, n22), nrow=2, ncol=2,
                  byrow=FALSE,
                  dimnames=list(
                    c("selected", "not"),
                    c(node, paste("NOT", node, sep=""))))
    ans <- .doHyperGInternal(numW, numB, numDrawn, numWdrawn, OVER)
    list(nodes=node,
         p.value=ans$p,
         method="Hypergeometric test using pyhper",
         observed=obs)
}

hg_test_factory <- function(selids, PCUT=0.05, COND=FALSE, OVER=TRUE)
{
    force(selids)
    function(start, g, prev_ans) {
        ## This is a global-level test using the Hypergeometric distribution.
        prev_ans_e <- new.env(hash=TRUE, parent=emptyenv())
        if (COND && length(prev_ans)) {
            nms <- sapply(prev_ans, function(x) x[["nodes"]])
            names(prev_ans) <- nms
            l2e(prev_ans, prev_ans_e)
        }
        
        lapply(unique(unlist(start)), function(aNode) {
            aNodeGenes <- geneIds(g, aNode)
            univ <- allGeneIds(g)
            if (COND) {
                ## XXX: this only makes sense with bottomup_iter for now.
                kids <- childrenOf(g, aNode)[[1]]
                kids_ans <- list()
                if (length(kids)) {
                    kids_ans <- mget(kids, prev_ans_e, ifnotfound=NA)
                    kids_ans <- kids_ans[sapply(kids_ans,
                                                function(x) !is.na(x[1]))]
                }
                if (length(kids_ans)) {
                    sigKids <- sapply(kids_ans, function(x) x[["nodes"]])
                    sigKids <- sigKids[sapply(kids_ans,
                                              function(x) x[["p"]] < PCUT)]
                    sigKidGenes <- unlist(lapply(sigKids,
                                                 function(x) geneIds(g, x)))
                    sigKidGenes <- unique(sigKidGenes)
                    aNodeGenes <- setdiff(aNodeGenes, sigKidGenes)
                    univ <- setdiff(univ, sigKidGenes)
                    selids <- setdiff(selids, sigKidGenes)
                }
            }
            numW <- length(aNodeGenes)
            numB <- length(univ) - numW
            numWdrawn <- sum(selids %in% aNodeGenes)
            numDrawn <- length(selids)

            make_hyperg_result(aNode, numW, numB, numDrawn,
                               numWdrawn, OVER)
        })
    }
}

## these are for testing the tree_iter functions
global_tdummy <- function(start, g, prev_ans)
{
    set.seed(0xeeff)
    lapply(unlist(start), function(n) {
        p <- runif(1)
        list(nodes=n, p=p)
    })
}

local_tdummy <- function(start, g, prev_ans)
{
    set.seed(0xeeff)
    lapply(start, function(theNodes) {
        p <- runif(1)
        list(nodes=theNodes, p=p)
    })
}

## basically, the test function tfun needs to accept a list.  For global
## mode, it will unlist and operate on individual nodes.  For local
## mode, it will operate on each element, assuming it is a "litter".

global_ndummy <- function(ans, g)
{
    ret <- lapply(ans, function(x) {
        if (x$p.value < 0.3)
          x$nodes
        else NULL
    })
    unlist(ret[!sapply(ret, is.null)])
}

ans_acceptor_factory <- function(PCUT=0.05)
{
    function(ans, g) {
        ret <- lapply(ans, function(x) {
            p <- x[["p.value"]]
            if (!is.na(p) && p < PCUT)
              x$nodes
            else NULL
        })
        unlist(ret[!sapply(ret, is.null)])
    }
}

simple_acceptor <- function(ans, g) {
    unlist(lapply(ans, function(x) x$nodes))
}

toDF <- function(ans) {
    lapply(ans, function(sub) {
        do.call(rbind, lapply(sub, function(r) {
            data.frame(bands=paste(r$nodes, collapse=", "),
                       p.value=r$p.value)
        }))
    })
}

cb_test <- function(selids, chrtree, level,
                    dir=c("up", "down"),
                    type=c("local", "global"),
                    pval=0.05,
                    conditional=FALSE)
{
    dir <- match.arg(dir)
    type <- match.arg(type)
    if (!(as.character(level) %in% names(chrtree@level2nodes)))
      stop("level must be one of ",
           paste(names(chrtree@level2nodes), collapse=", "))
    if ((type == "local" || dir == "down") && conditional)
      stop('conditional can only be used for type="global" and dir="up"')
    if (pval < 0 || pval > 1)
      stop("'pval' must be bewteen 0 and 1")

    iter <- switch(dir,
                   up=bottomup_iter,
                   down=topdown_iter,
                   stop("'dir' must be 'up' or 'down'"))
    tfun <- switch(type,
                   local=local_test_factory(selids, pval),
                   global=hg_test_factory(selids, pval, COND=conditional),
                   stop("'type' must be 'local' or 'global'"))

    start <- switch(type,
                    local={
                        oneUp <- as.character(level -1L)
                        s <- childrenOf(chrtree, chrtree@level2nodes[[oneUp]])
                        s[listLen(s) > 0]
                    },
                    global=chrtree@level2nodes[[level]])

    iter(chrtree, start, tfun, ans_acceptor_factory(pval))
}

## notes: return the table because then odds ratio becomes a method,
## e.g. as expected.  This also starts to unify results from global and
## local tests.  oddratio will give NA if k x 2 table.
