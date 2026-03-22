#' Write the 10X molecule information file
#'
#' Write a data frame of molecule information into a HDF5 file.
#'
#' @param path String containing the path to the file.
#' @param data List returned by \code{\link{readMolecules}} or \code{\link{simulateMolecules}}.
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
#' @return For \code{writeMolecules}, a new HDF5 file is created at \code{path} containing the molecule information.
#' \code{NULL} is invisibly returned.
#'
#' For \code{simulateMolecules}, a list is returned containing \code{molecules}, \code{features}, and (for \code{version="3"}) \code{libraries}.
#' Each of these mimic the corresponding entry in the list returned by \code{\link{readMolecules}}.
#'
#' @details
#' These functions are primarily intended for examples and testing of \code{readMolecules}.
#' The files written by \code{writeMolecules} may not be sufficiently compliant with 10X's specification for use in external applications.
#' If you are running into interoperability issues, file an issue at \url{https://github.com/LTLA/crio} and we'll see what we can do.
#' 
#' @author Aaron Lun
#'
#' @examples
#' sim <- simulateMolecules()
#' tmp <- tempfile(fileext=".h5")
#' writeMolecules(tmp, sim)
#' rhdf5::h5ls(tmp)
#' 
#' @export
#' @importFrom rhdf5 h5write h5createGroup h5createFile
#' @importFrom jsonlite toJSON
writeMolecules <- function(
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

    stopifnot(!is.null(molecules$barcode))
    if (version == "2") {
        barcodes <- as.character(molecules$barcode)
        barcodes <- sub("-.*", "", barcodes) # remove any gem groups.
        write_cell_barcodes(barcodes, temp, "barcode")
    } else {
        actual.barcodes <- as.factor(molecules$barcode) # convert it to a factor but only if it wasn't already one.
        h5write_type(as.integer(actual.barcodes) - 1L, temp, group = NULL, name = "barcode_idx", type = "H5T_NATIVE_UINT64") 
        h5write(levels(actual.barcodes), temp, "barcodes")
    }

    stopifnot(!is.null(molecules$feature))
    if (version == "2") {
        feature.field <- "gene"
    } else {
        feature.field <- "feature_idx"
    }
    h5write_type(molecules$feature - 1L, temp, group = NULL, name = feature.field, type = "H5T_NATIVE_UINT32")

    stopifnot(!is.null(molecules$count))
    if (version == "2") {
        read.field <- "reads"
    } else {
        read.field <- "count"
    }
    h5write_type(molecules$count, temp, group = NULL, name = read.field, type = "H5T_NATIVE_UINT32")

    stopifnot(!is.null(molecules$umi))
    h5write_type(molecules$umi, temp, group = NULL, name = "umi", type = "H5T_NATIVE_UINT32")
    stopifnot(!is.null(molecules$gem.group))
    h5write_type(molecules$gem.group, temp, group = NULL, name = "gem_group", type = "H5T_NATIVE_UINT16")

    stopifnot(!is.null(features$id))
    if (version == "2") {
        h5write(features$id, temp, "gene_ids")
    } else {
        stopifnot(!is.null(features$type))
        h5createGroup(temp, "features")
        h5write(features$id, temp, "features/id")
        h5write(features$type, temp, "features/feature_type")
    }

    if (version == "3") {
        stopifnot(!is.null(molecules$library))
        h5write_type(molecules$library - 1L, temp, group = NULL, name = "library_idx", type = "H5T_NATIVE_UINT16")
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
#' @rdname writeMolecules
simulateMolecules <- function(
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

    mol.info <- DataFrame(barcode = factor(barcodes), umi = umi, gem.group = gem_group, feature = feature, count = reads)
    feature.info <- DataFrame(id = sprintf("FEATURE-%s", seq_len(num.features)))
    output <- list(molecules = mol.info, features = feature.info)

    if (version == "3") {
        output$features$type <- rep("Gene Expression", nrow(output$features))
        output$molecules$library <- mol.info$gem.group
        output$libraries <- lapply(seq_len(num.gem.groups), function(i) {
            list(library_type = "unknown", library_id = i - 1L, gem_group = i)
        })
    }

    output

}
