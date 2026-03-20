# library(testthat); library(crio); source("test-encodeSequences.R")

test_that("encodeSequences works as expected", {
    truth <- c("ACGT", "AATT", "GGCC", "TGCA")
    enc <- encodeSequences(truth)
    expect_identical(truth, decodeSequences(enc, lengths=4))

    expect_error(encodeSequences(strrep("A", 50)), "15 nt")
    expect_error(encodeSequences("ABCDE"), "unrecognized")

    test <- decodeSequences(enc, lengths=c(3,4,5,6))
    expect_identical(test, c("CGT", "AATT", "AGGCC", "AATGCA"))
})
