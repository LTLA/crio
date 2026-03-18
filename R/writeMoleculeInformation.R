#' Write the 10X molecule information file
#'
#' Write a data frame of molecule information into a HDF5 file.
#'
#' @param path String containing the path to the file.
#' @param data List returned by \code{\link{readMoleculeInformation}} or \code{\link{simulateMoleculeInformation}}.
#' This should contain the \code{"molecules"}, \code{"features"} and (for \code{version="3"}) \code{"libraries"} entries.
#' @param version String specifying the version of the 10X molecule information format to write.
#' This can be either \code{"2"} or \code{"3"}.
#' @param num.molecules Integer specifying the number of molecules to simulate.
#' @param num.features Integer specifying the number of features in the simulation.
#' @param barcode.length Integer specifying the length of the cell barcodes.
#' @param umi.length Integer specifying the length of the UMI.
#' @param ave.num.reads Number specifying the average of the number of reads per molecule.
#' @param num.gem.groups Integer specifying the number of GEM groups. 
#' This is also treated as the number of libraries for \code{version="3"}.
#'
#' @return For \code{writeMoleculeInformation}, a new HDF5 file is created at \code{path} containing the molecule information.
#' \code{NULL} is invisibly returned.
#'
#' For \code{simulateMoleculeInformation}, a list is returned containing \code{molecules}, \code{features}, and (for \code{version="3"}) \code{libraries}.
#' Each of these mimic the corresponding entry in the list returned by \code{\link{readMoleculeInformation}}.
#'
#' @details
#' These functions are primarily intended for examples and testing of \code{readMoleculeInformation}.
#' The files written by \code{writeMoleculeInformation} not be sufficiently compliant with 10X's specification for use in external applications.
#' If you are running into interoperability issues, file an issue at \url{https://github.com/LTLA/crio} and we'll see what we can do.
#' 
#' @author Aaron Lun
#'
#' @examples
#' sim <- simulateMoleculeInformation()
#' tmp <- tempfile(fileext=".h5")
#' writeMoleculeInformation(tmp, sim)
#' rhdf5::h5ls(tmp)
#' 
#' @export
#' @importFrom rhdf5 h5write h5createGroup h5createFile
#' @importFrom jsonlite toJSON
writeMoleculeInformation <- function(
    path,
    data,
    version = c("3", "2")
) {
    temp <- tempfile(tmpdir = dirname(path))
    on.exit({
        if (file.exists(temp)) {
            unlink(temp)
        }
    })
    h5 <- h5createFile(temp)

    version <- match.arg(version)
    molecules <- data$molecules
    features <- data$features
    libraries <- data$libraries

    stopifnot(!is.null(molecules$Barcode))
    if (version == "2") {
        barcodes <- as.character(molecules$Barcode)
        barcodes <- sub("-.*", "", barcodes) # remove any gem groups.
        write_cell_barcodes(barcodes, temp, "barcode")
    } else {
        actual.barcodes <- as.factor(molecules$Barcode) # convert it to a factor but only if it wasn't already one.
        h5write_type(as.integer(actual.barcodes) - 1L, temp, group = NULL, name = "barcode_idx", type = "H5T_NATIVE_UINT64") 
        h5write(levels(actual.barcodes), temp, "barcodes")
    }

    stopifnot(!is.null(molecules$Feature))
    if (version == "2") {
        feature.field <- "gene"
    } else {
        feature.field <- "feature_idx"
    }
    h5write_type(molecules$Feature - 1L, temp, group = NULL, name = feature.field, type = "H5T_NATIVE_UINT32")

    stopifnot(!is.null(molecules$Count))
    if (version == "2") {
        read.field <- "reads"
    } else {
        read.field <- "count"
    }
    h5write_type(molecules$Count, temp, group = NULL, name = read.field, type = "H5T_NATIVE_UINT32")

    stopifnot(!is.null(molecules$Umi))
    h5write_type(molecules$Umi, temp, group = NULL, name = "umi", type = "H5T_NATIVE_UINT32")
    stopifnot(!is.null(molecules$GemGroup))
    h5write_type(molecules$GemGroup, temp, group = NULL, name = "gem_group", type = "H5T_NATIVE_UINT16")

    stopifnot(!is.null(features$ID))
    if (version == "2") {
        h5write(features$ID, temp, "gene_ids")
    } else {
        stopifnot(!is.null(features$Type))
        h5createGroup(temp, "features")
        h5write(features$ID, temp, "features/id")
        h5write(features$Type, temp, "features/feature_type")
    }

    if (version == "3") {
        stopifnot(!is.null(molecules$Library))
        h5write_type(molecules$Library - 1L, temp, group = NULL, name = "library_idx", type = "H5T_NATIVE_UINT16")
        h5write(as.character(toJSON(libraries, auto_unbox=TRUE)), temp, "library_info")
    }

    # Only moving at the end to avoid messing up 'path' if any of the previous steps were unsuccessful.
    unlink(path)
    file.rename(temp, path)
    invisible(NULL)
}

#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom stats rpois
#' @rdname writeMoleculeInformation
simulateMoleculeInformation <- function(
    num.molecules = 10000,
    barcode.length = 4, 
    num.features = 20,
    num.gem.groups = 1,
    umi.length = 10,
    ave.num.reads = 10,
    version = c("3", "2")
) {
    version <- match.arg(version)

    barcodes <- ""
    for (y in seq_len(barcode.length)) {
        barcodes <- paste0(barcodes, sample(c("A", "C", "G", "T"), num.molecules, replace=TRUE)) 
    }

    umi <- sample(4L ^ as.integer(umi.length), num.molecules) - 1L
    feature <- sample(num.features + 1L, num.molecules, replace = TRUE)
    reads <- pmax(1L, rpois(num.molecules, lambda = ave.num.reads))

    gem_group <- sample(num.gem.groups, num.molecules, replace=TRUE)
    if (version == "3") {
        barcodes <- sprintf("%s-%s", barcodes, gem_group)
    }

    mol.info <- DataFrame(Barcode = factor(barcodes), Umi = umi, GemGroup = gem_group, Feature = feature, Count = reads)
    feature.info <- DataFrame(ID = sprintf("FEATURE-%s", seq_len(num.features)))
    output <- list(molecules = mol.info, features = feature.info)

    if (version == "3") {
        output$features$Type <- rep("Gene Expression", nrow(output$features))
        output$molecules$Library <- mol.info$GemGroup
        output$libraries <- lapply(seq_len(num.gem.groups), function(i) {
            list(library_type = "unknown", library_id = i - 1L, gem_group = i)
        })
    }

    output

}
