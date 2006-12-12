library("hgu95av2")
quiet <- TRUE

testValid <- function() {
    univ <- 1:100
    sel <- sample(univ, 10)
    p1 <- new("KEGGHyperGParams",
                       geneIds=sel,
                       universeGeneIds=univ,
                       annotation="hgu95av2",
                       pvalueCutoff=0.05,
                       testDirection="over")

    p2 <- new("GOHyperGParams",
                       geneIds=sel,
                       universeGeneIds=univ,
                       annotation="hgu95av2",
                       pvalueCutoff=0.05,
                       testDirection="over")
}

testDupGeneIds <- function() {
    univ <- 1:100
    sel <- as.integer(c(1, 1, 1, 2, 3, 4, 5))
    checkException(new("KEGGHyperGParams",
                       geneIds=sel,
                       universeGeneIds=univ,
                       annotation="hgu95av2",
                       pvalueCutoff=0.05,
                       testDirection="over"),
                   silent=quiet)

    checkException(new("GOHyperGParams",
                       geneIds=sel,
                       universeGeneIds=univ,
                       annotation="hgu95av2",
                       pvalueCutoff=0.05,
                       testDirection="over"),
                   silent=quiet)

}

testDupUniverseGeneIds <- function() {
    univ <- rep(1:100, 2)
    sel <- 1:10
    checkException(new("KEGGHyperGParams",
                       geneIds=sel,
                       universeGeneIds=univ,
                       annotation="hgu95av2",
                       pvalueCutoff=0.05,
                       testDirection="over"),
                   silent=quiet)

    checkException(new("GOHyperGParams",
                       geneIds=sel,
                       universeGeneIds=univ,
                       annotation="hgu95av2",
                       pvalueCutoff=0.05,
                       testDirection="over"),
                   silent=quiet)
}

testGeneIdsMatchUniverseGeneIds <- function() {
    univ <- 1:100
    sel <- as.character(1:10)
    checkException(new("KEGGHyperGParams",
                       geneIds=sel,
                       universeGeneIds=univ,
                       annotation="hgu95av2",
                       pvalueCutoff=0.05,
                       testDirection="over"),
                   silent=quiet)

    checkException(new("GOHyperGParams",
                       geneIds=sel,
                       universeGeneIds=univ,
                       annotation="hgu95av2",
                       pvalueCutoff=0.05,
                       testDirection="over"),
                   silent=quiet)
}

testBadPvalueCutoff <- function() {
    univ <- as.character(1:100)
    sel <- as.character(1:10)
    checkException(new("KEGGHyperGParams",
                       geneIds=sel,
                       universeGeneIds=univ,
                       annotation="hgu95av2",
                       pvalueCutoff=5,
                       testDirection="over"),
                   silent=quiet)

    checkException(new("GOHyperGParams",
                       geneIds=sel,
                       universeGeneIds=univ,
                       annotation="hgu95av2",
                       pvalueCutoff=-.05,
                       testDirection="over"),
                   silent=quiet)
}
