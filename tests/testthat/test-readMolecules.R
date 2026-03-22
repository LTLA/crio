# library(crio); library(testthat); source("test-readMolecules.R")

set.seed(909)
test_that("readMolecules works correctly", {
    tmp <- tempfile(fileext=".h5")
    simulated <- simulateMolecules()
    writeMolecules(tmp, simulated)
    roundtrip <- readMolecules(tmp) 

    keep <- simulated$molecules$feature <= nrow(simulated$features)
    expect_identical(simulated$molecules[keep,], roundtrip$molecules)
    expect_identical(simulated$features, roundtrip$features)
    expect_identical(simulated$libraries, roundtrip$libraries)

    roundtrip <- readMolecules(tmp, keep.unmapped=TRUE) 
    expect_identical(simulated$molecules, roundtrip$molecules)
})

set.seed(9091)
test_that("partial extraction in readMolecules works correctly", {
    tmp <- tempfile(fileext=".h5")
    simulated <- simulateMolecules()
    writeMolecules(tmp, simulated)

    keep <- simulated$molecules$feature <= nrow(simulated$features)
    full <- readMolecules(tmp)

    # Discounting the GEM.
    subbed <- readMolecules(tmp, get.gem=FALSE)
    copy <- simulated
    copy$molecules <- copy$molecules[keep,]
    copy$molecules$gem.group <- NULL
    expect_identical(subbed, copy)

    # Discounting the genes. 
    subbed <- readMolecules(tmp, get.feature=FALSE)
    copy <- simulated
    copy$molecules <- copy$molecules[keep,]
    copy$molecules$feature <- NULL
    expect_identical(subbed, copy)

    fullun <- readMolecules(tmp, keep.unmapped=TRUE)
    subbed <- readMolecules(tmp, get.feature=FALSE, keep.unmapped=TRUE)
    fullun$molecules$feature <- NULL
    expect_identical(subbed, fullun)

    # Discounting everything.
    subbed <- readMolecules(tmp, get.gem=FALSE, get.count=FALSE, get.barcode=FALSE, get.feature=FALSE, get.umi=FALSE, get.library=FALSE)
    expect_identical(nrow(subbed$molecules), nrow(full$molecules))
    expect_identical(ncol(subbed$molecules), 0L)
    expect_identical(subbed$features, full$features)

    # Discounting everyting but keeping unmapped reads.
    subbed <- readMolecules(tmp, get.gem=FALSE, get.count=FALSE, get.barcode=FALSE, get.feature=FALSE, get.umi=FALSE, get.library=FALSE, keep.unmapped=TRUE)
    expect_identical(nrow(subbed$molecules), nrow(fullun$molecules))
    expect_identical(ncol(subbed$molecules), 0L)
    expect_identical(subbed$features, full$features)

    # Skipping the library information.
    nolib <- readMolecules(tmp, extract.library.info = FALSE)
    expect_null(nolib$libraries)
})

set.seed(908)
test_that("readMolecules works with older CellRanger versions", {
    simulated <- simulateMolecules(version="2")
    expect_null(simulated$libraries)

    tmp <- tempfile(fileext=".h5")
    writeMolecules(tmp, simulated, version="2")
    restored <- readMolecules(tmp, keep.unmapped=TRUE) # auto-detects the version.
    expect_identical(restored, simulated)

    restored <- readMolecules(tmp, keep.unmapped=TRUE, barcode.length=4) # also works if we pass the barcode length.
    expect_identical(restored, simulated)

    expect_error(readMolecules(tmp, keep.unmapped=TRUE, barcode.length=-1), 'non-negative')
    expect_error(readMolecules(tmp, keep.unmapped=TRUE, barcode.length=c(1,2,3)), 'scalar')

    lowered <- restored
    lowered$molecules$barcode <- tolower(lowered$molecules$barcode)
    restored <- readMolecules(tmp, keep.unmapped=TRUE)
    expect_identical(restored, simulated)

    # If we accidentally give it version 3 and write it as version 2, it's still fine.
    set.seed(9999)
    simulated <- simulateMolecules()

    tmp <- tempfile(fileext=".h5")
    writeMolecules(tmp, simulated, version="2")
    restored <- readMolecules(tmp, keep.unmapped=TRUE) # auto-detects the version.
    restored$molecules$barcode <- factor(paste0(restored$molecules$barcode, "-1"))

    copy <- simulated 
    copy$molecules$library <- NULL
    copy$features$type <- NULL
    copy$libraries <- NULL
    expect_identical(restored, copy)

    # Checking various error conditions. 
    failed <- simulated
    levels(failed$molecules$barcode)[1] <- ""
    ftmp <- tempfile()
    expect_error(writeMolecules(ftmp, failed, version="2"))

    failed <- simulated
    levels(failed$molecules$barcode)[1] <- strrep("A", 50)
    ftmp <- tempfile()
    expect_error(writeMolecules(ftmp, failed, version="2"))

    failed <- simulated
    levels(failed$molecules$barcode)[1] <- "whee"
    ftmp <- tempfile()
    expect_error(writeMolecules(ftmp, failed, version="2"))
})

set.seed(907)
test_that("readMolecules works with silly inputs containing no molecules", {
    simulated <- simulateMolecules(num.molecules=0)
    temp <- tempfile(fileext=".h5")
    writeMolecules(temp, simulated)

    roundtrip <- readMolecules(temp)
    expect_identical(nrow(roundtrip$molecules), 0L)
    expect_gt(nrow(roundtrip$features), 0)

    # Checking that it doesn't throw up with automatic barcode detection.
    temp2 <- tempfile(fileext=".h5")
    writeMolecules(temp2, simulated, version="2")
    roundtrip2 <- readMolecules(temp2, version="2")
    expect_identical(nrow(roundtrip2$molecules), 0L)

    # Checking that it behaves when there aren't even any genes. 
    simulated0 <- simulateMolecules(num.molecules=0, num.features=0)
    temp0 <- tempfile(fileext=".h5")
    writeMolecules(temp0, simulated0, version="2")
    roundtrip0 <- readMolecules(temp0)
    expect_identical(nrow(roundtrip0$molecules), 0L)
    expect_identical(nrow(roundtrip0$features), 0L)
})
