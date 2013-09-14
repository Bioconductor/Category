setMethod(hyperg, "list",
    function(assayed, significant, universe, representation=c("over", "under"),
             ...)
{
    .hyperg <- selectMethod(hyperg, "character")
    result <- Map(.hyperg, assayed, significant,
        MoreArgs=list(universe=universe, representation=representation, ...))
    as.data.frame(t(simplify2array(result)))
})

setMethod(hyperg, "character",
     function(assayed, significant, universe, representation=c("over", "under"))
{
    if (!all(assayed %in% universe))
        stop("some 'assayed' genes not in 'universe'")
    if (!all(significant %in% assayed))
        stop("some 'significant' genes not in 'assayed'")
    representation <- match.arg(representation)
    test <- representation== "over"

    white.balls.drawn <- length(intersect(significant, assayed))
    white.balls.in.urn <- length(assayed)
    total.balls.in.urn <- length(universe)
    black.balls.in.urn <- total.balls.in.urn - white.balls.in.urn
    balls.pulled.from.urn <- length(assayed)

    c(assayed=length(assayed), significant=length(significant),
      universe=length(universe), representation=representation,
      .doHyperGInternal(white.balls.in.urn, black.balls.in.urn,
                        balls.pulled.from.urn, white.balls.drawn,
                        test))
})
