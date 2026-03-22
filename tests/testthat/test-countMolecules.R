# library(crio); library(testthat); source("test-countMolecules.R")

library(Matrix)
REFERENCE <- function(molinfo, num.features, use.count=FALSE) {
    if (use.count) {
        vals <- molinfo$count
    } else {
        vals <- rep(1, nrow(molinfo))
    }
    keep <- molinfo$feature <= num.features

    sparseMatrix(
        i = molinfo$feature[keep],
        j = as.integer(molinfo$barcode)[keep],
        x = vals[keep],
        dims = c(num.features, nlevels(molinfo$barcode)),
        dimnames = list(NULL, levels(molinfo$barcode))
    )
}

set.seed(123)
test_that("countMolecules works correctly", {
    sim <- simulateMoleculeInformation(num.features=20)
    expect_true(any(sim$molecules$feature == 21)) # check that there were, indeed, features beyond the range.

    out <- countMolecules(20, sim$molecules$feature, sim$molecules$barcode)
    ref <- REFERENCE(sim$molecules, num.features=20)
    expect_equal(ref, out)

    # Same result if we removed the offending entries manually.
    keep <- sim$molecules$feature <= 20
    out.filtered <- countMolecules(20, sim$molecules$feature[keep], sim$molecules$barcode[keep])
    expect_equal(ref, out.filtered)

    # Works with a SVT_SparseMatrix.
    out.svt <- countMolecules(20, sim$molecules$feature, sim$molecules$barcode, class="SVT_SparseMatrix")
    ref.svt <- as(ref, "SVT_SparseArray")
    type(ref.svt) <- "integer"
    expect_identical(out.svt, ref.svt)

    # Works with multiple threads.
    out.parallel <- countMolecules(20, sim$molecules$feature, sim$molecules$barcode, num.threads=2)
    expect_equal(ref, out.parallel)
    out.svt.parallel <- countMolecules(20, sim$molecules$feature, sim$molecules$barcode, class="SVT_SparseMatrix", num.threads=2)
    expect_equal(ref.svt, out.svt.parallel)
})

set.seed(456)
test_that("countMolecules filters out NA features and barcodes", {
    sim <- simulateMoleculeInformation(num.features=20)
    sim$molecules$feature[1] <- NA
    sim$molecules$barcode[2] <- NA
    out <- countMolecules(20, sim$molecules$feature, sim$molecules$barcode)

    remaining <- 3:nrow(sim$molecules)
    ref <- countMolecules(20, sim$molecules$feature[remaining], sim$molecules$barcode[remaining])
    expect_equal(ref, out)
})

set.seed(789)
test_that("countMolecules sums the counts", {
    sim <- simulateMoleculeInformation(num.features=20)

    out <- countMolecules(20, sim$molecules$feature, sim$molecules$barcode, count=sim$molecules$count)
    ref <- REFERENCE(sim$molecules, num.features=20, use.count=TRUE)
    expect_equal(ref, out)

    # Works with a SVT_SparseMatrix.
    out.svt <- countMolecules(20, sim$molecules$feature, sim$molecules$barcode, class="SVT_SparseMatrix", count=sim$molecules$count)
    ref.svt <- as(ref, "SVT_SparseArray")
    type(ref.svt) <- "integer"
    expect_identical(out.svt, ref.svt)
})
