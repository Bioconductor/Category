qRequire <- function(pkg) {
    suppressWarnings(require(pkg, character.only=TRUE,
                             quietly=TRUE, warn.conflicts=FALSE))
}

pvalFromPermMat <- function(obs, perms) {
    N <- ncol(perms)
    pvals <- matrix(as.double(NA), nr=nrow(perms), ncol=2)
    dimnames(pvals) <- list(rownames(perms), c("Lower", "Upper"))

    tempObs <- rep(obs, ncol(perms))
    dim(tempObs) <- dim(perms)
    pvals[ , 1] <- rowSums(perms <= tempObs)/N
    pvals[ , 2] <- rowSums(perms >= tempObs)/N
    pvals
}

gseaperm <- function(eset, fac, mat, nperm) {
    mkSparseMat <-
      if (qRequire("Matrix")) function(x) Matrix(x, sparse=TRUE)
      else matrix

    geneNames <- colnames(mat)
    if (is.null(geneNames))
      stop("'mat' argument must have column names")
    eset <- eset[colnames(mat), ]
    if (nrow(eset) != ncol(mat))
      warning("'eset' and 'mat' genes not identical")
    if (nrow(eset) < 2)
      stop("need at two genes in common between 'eset' and 'mat'")
    cAmat <- mkSparseMat(mat)

    obs <- rowttests(eset, fac, tstatOnly=TRUE)[["statistic"]]
    obs <- as.vector(cAmat %*% obs)

    permMat <- matrix(0, nrow=nrow(eset), ncol=nperm)
    i <- 1L
    while (i < (nperm + 1)) {
        p1 <- sample(fac)
        permMat[ , i] <- rowttests(eset, p1, tstatOnly=TRUE)[["statistic"]]
        i <- i + 1L
    }
    permMat <- cAmat %*% permMat

    pvalFromPermMat(obs, permMat)
}
