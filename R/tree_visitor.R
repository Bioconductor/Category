tree_iter <- function(g, start, tfun, nfun, relationOf)
{
    i <- 1L
    e <- new.env(hash=TRUE, parent=emptyenv())
    prev_ans <- list()
    while (TRUE) {
        if (!length(start)) break
        ans <- tfun(start, g, prev_ans)
        k <- as.character(i)
        e[[k]] <- ans
        prev_ans <- ans
        i <- i + 1L
        accepted <- unique(nfun(ans, g))
        if (length(accepted))
          start <- relationOf(g, accepted)
        else
          start <- character(0)
    }
    as.list(e)
}

topdown_iter <- function(g, start, tfun, nfun)
{
    tree_iter(g, start, tfun, nfun, childrenOf)
}

bottomup_iter <- function(g, start, tfun, nfun)
{
    tree_iter(g, start, tfun, nfun, parentOf)
}
