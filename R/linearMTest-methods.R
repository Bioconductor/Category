
### no work done on this yet.  I think I'll start by writing prototype
### code that works for the Cytoband paper example, and then adapt
### that.

setMethod("linearMTest",
          signature(p = "LinearMParams"),
          function(p) .linearMTestInternal(p))



.linearMTestInternal <- function(p, className="LinearMResult")
{

    ## in the hyperG case, this function is used as a wrapper that
    ## calls different backends based on class (sort of a fake method
    ## dispatch).  As of now, we have LinearM code only for the
    ## chromosome band case, so we won't pay too much attention to
    ## making the code optimised for future extensions.

    ## we don't use categoryToEntrezBuilder(), because our code so far
    ## is only meant for nested categories.


    if (categoryName(p) != "ChrMap")
        stop("'linearMTest' currently only works for chromosome band categories")
    className <- paste(categoryName(p), className, sep = "")
    ## p <- makeValidParams(p) # FIXME: not written yet

    if (numNodes(p@chrGraph) == 0)
        p@chrGraph <- makeChrBandGraph(p@annotation, p@universeGeneIds)
    ## Note: unlike the hyperG case, we do not store results of the
    ## tests in a graph. This may be a good idea, but needs thought
    ## because we potentially remove some nodes from consideration if
    ## they have too few genes.
    ans <-
        linearMTest_ChrMap(stats = p@geneStats,
                           g = p@chrGraph,
                           upreg = p@testDirection == "up",
                           chip = annotation(p),
                           cutoff = pvalueCutoff(p),
                           global = TRUE,
                           min.genes = p@minSize,
                           conditional = conditional(p))
    new(className,
        pvalues = if (p@conditional) ans$conditional.pvals else ans$marginal.pvals,
        annotation = p@annotation,
        geneIds = p@universeGeneIds,
        testName = p@categoryName,
        pvalueCutoff = p@pvalueCutoff,
        testDirection = p@testDirection,
        minSize = p@minSize,
        pvalue.order = order(if (p@conditional) ans$conditional.pvals else ans$marginal.pvals),
        conditional = p@conditional,
        chrGraph = p@chrGraph,
        catToGeneId = ans$catToGeneId)
}



## Here's the background for determining conditional significance of
## chromosome bands: Assume that we know how to determine marginal
## significance (this is done by fitting the model t ~ 1 + x, where t
## is the vector of t-stats and x is the indicator of category
## membership, where rows are genes or probesets).  Let's say we have
## a graph that represents the chromosome band structure, and an
## incidence matrix giving membership of genes in chromosome bands.
## Start with graph induced by marginally (as opposed to
## conditionally) significant categories.  Alternatively, start with
## full graph.  Declare significant leaves (there cannot be
## non-significant leaves in the first case, I think) to be
## conditionally significant.  For every other marginally significant
## node, assume conditional significance has been determined for all
## children (if not, do those first).  Fit a linear model with t-stats
## as response, and indicators for all CONDITIONALLY (not marginally)
## significant children as predictors, as well as an indicator for the
## current category.  If the current category is significant in the
## result, declare current node as conditionally significant.

## This approach is also valid for non-nested genesets such as GO, but
## then we will need to worry about overlap between rest (which we
## could try to do pairwise).  This is not an issue with chromosome
## bands.

## ok, so how to do this?

## First, we need a method for finding all children of a node in a
## graph.  We'll use a quick and dirty solution that repeatedly
## follows edges until there are no more to find

.getChildrenInGraph <- function(g, init)
    ## g: graph
    ## init: vector (typically length-1?) of node names
{
    gm <- as(g, "matrix")
    final <- logical(nrow(gm))
    names(final) <- rownames(gm)
    final[init] <- TRUE
    done <- FALSE
    while (!done) {
        ncur <- sum(final)
        direct.children <-
            apply(gm[final, , drop = FALSE], 1,
                  function(x) which(x != 0))
        if (is.list(direct.children))
            direct.children <- unlist(direct.children)
        final[direct.children] <- TRUE
        nnew <- sum(final)
        done <- nnew == ncur
    }
    names(final)[final]
}



## With this tool, we're now ready to do the sequence of tests.
## Basically, we can test a node if it has no children.  If it does,
## then the children must have been tested first, and those that are
## significant would be included in the model.

## This function does the (trivial) marginal test for one category

linearMTestMarginal <- function(stats, x, upreg = TRUE)
{
    lm.res <- lm(stats ~ 1 + x)
    tt <- summary(lm.res)$coef[2, 3]
    if (upreg) 1 - pnorm(tt) else pnorm(tt)
}


## This function does the conditional test (for the whole graph, not
## just one category).  As the marginal/unconditional tests are done
## anyway, we allow skipping the conditional testing through an
## argument.

linearMTest_ChrMap <-
    function(stats, g, upreg = TRUE, chip = "hgu95av2",
             cutoff = 0.01,
             global = TRUE,
             min.genes = 4,
             conditional = TRUE,
             verbose = getOption("verbose"))
    ## upreg: whether looking for H_1: mu > 0 or mu < 0
    ## global: global means all genes in 'stats' considered universe.
    ##    local means only ones which belong in some node in g.  
    ##    FIXME: not sure if this is the correct interpretation.
    ## min.genes: keep only nodes with these many genes
{
    ## print(system.time(
    inc.mat <- t(MAPAmat(chip, univ = names(stats)))
    ## ))
    stopifnot(setequal(names(stats), rownames(inc.mat)))

    ## leave out nodes that have no genes (or rather, too few genes)
    cb.names <- colnames(inc.mat)
    cb.names <- cb.names[ (colSums(inc.mat) >= min.genes) & (cb.names %in% nodes(g)) ]
    ## cb.names <- cb.names[grep("^[1-9XY]", cb.names)]
    g <- subGraph(cb.names, g)
    inc.mat <- inc.mat[, cb.names]
    gene.names <- rownames(inc.mat)
    stats <- stats[gene.names]
    catToGeneId <- vector(mode = "list", length = length(cb.names))
    names(catToGeneId) <- cb.names
    for (nm in cb.names)
        catToGeneId[[nm]] <- gene.names[which(inc.mat[, nm] == 1)]

    if (!global)
    {
        keep <- rowSums(inc.mat) > 0
        stats <- stats[keep]
        inc.mat <- inc.mat[keep, ]
    }

    ## set up logical vectors to track quantities of interest
    marginal.pvals <- ## implies names(.) <- cb.names
        sapply(cb.names,
               function(n) {
                   linearMTestMarginal(stats, inc.mat[, n],
                                       upreg = upreg)
               })
    marginal <- marginal.pvals < cutoff
    if (!conditional)
        return(list(marginal = marginal,
                    marginal.pvals = marginal.pvals,
                    catToGeneId = catToGeneId))

    ## otherwise, proceed with conditional testing
    excluded <- conditional <- logical(length(cb.names))
    names(excluded) <- names(conditional) <- cb.names

    ## containter for adjusted (conditional) p-values
    adjusted.pvals <- marginal.pvals

    ## first iteration, figure out which ones have no children, i.e.,
    ## the ones that are leaves
    gl <- leaves(g, "out")
    conditional[gl] <- adjusted.pvals[gl] < cutoff
    excluded[gl] <- TRUE

    while (!all(excluded))
    {
        subg <- subGraph(cb.names[!excluded], g)
        subl <- leaves(subg, "out")
        if (verbose) cat(subl, "\n")

        ## loop over leaves and decide if they are conditionally
        ## significant

        for (i in subl)
        {
            ## FIXME: does it make sense to check conditional
            ## significance when not marginally significant?
            if (!is.na(marginal[i]) && marginal[i]) ## else nothing to do
            {
                if (verbose) cat("Processing node: ", i, "\n")
                ## need to find children that are conditionally
                ## significant

                ## how to find children?  Use function .getChildrenInGraph.
                ## Work on a smaller graph: current + excluded.  This
                ## is not necessary, but is done in the hope that this
                ## is more efficient.  FIXME: Don't know that for a
                ## fact, should profile.
                smallg <- subGraph(c(i, cb.names[excluded]), g)
                ## get all children of i.  This always includes i
                all.children <- .getChildrenInGraph(smallg, i)
                ## keep only those that are signif. This step
                ## automtically excludes i since conditional[i] is
                ## guaranteed to be FALSE at this point
                sig.children <- all.children[conditional[all.children]]

                ## Great. Now to fit a linear model.  We already have
                ## the response.  The design matrix should have (1) an
                ## intercept (unless we think the response should be
                ## centered) (2) one column for each sig.children and
                ## (3) one column for the current node

                ## this is a bit inefficient, but easiest on the mind
                ## (the more efficient alternative would be to use
                ## lm.fit)
                design.mat <-
                    cbind(y = stats,
                          inc.mat[, c(sig.children, i), drop = FALSE])
                design.mat <- as.data.frame(design.mat)
                lm.res <- lm(formula(design.mat), design.mat)
                cond.t <- (summary(lm.res)$coef[, 3])[ncol(design.mat)]
                cond.pval <- if (upreg) 1 - pnorm(cond.t) else pnorm(cond.t)

                ## good, now set these results
                adjusted.pvals[i] <- cond.pval
                conditional[i] <- !is.na(cond.pval) && cond.pval < cutoff
            }
        }
        ## done with all the leaves, so exclude them (and move on to
        ## next layer)
        excluded[subl] <- TRUE
    }
    list(marginal = marginal,
         conditional = conditional,
         marginal.pvals = marginal.pvals,
         conditional.pvals = adjusted.pvals,
         catToGeneId = catToGeneId)
}


