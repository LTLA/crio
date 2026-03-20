#' Read the 10X molecule information file
#' 
#' Extract relevant fields from the molecule information HDF5 file.
#' 
#' @param path String containing the path to the molecule information HDF5 file.
#' @param barcode.length Integer specifying the length of the cell barcode.
#' If \code{NULL}, the barcode length is guessed from the 2-bit encoded integers. 
#' Only relevant when \code{version="2"}.
#' @param keep.unmapped Boolean indicating whether unmapped molecules should be reported.
#' @param get.barcode,get.umi,get.gem,get.feature,get.count,get.library 
#' Logical scalar indicating whether the corresponding field should be extracted for each molecule.
#' @param version String specifying the version of the 10X molecule information format to read data from.
#' This can be either \code{"2"} or \code{"3"}.
#' If \code{NULL}, the version is auto-detected from the contents of \code{path}.
#' @param extract.library.info Boolean whether the library information should be extracted.
#' Only relevant when \code{version="3"}.
#' 
#' @return A named list containing molecule information.
#'
#' The first entry of this list is \code{molecules}, a \link[S4Vectors]{DataFrame} where each row corresponds to a single transcript molecule.
#' This contains the following fields:
#' \describe{
#' \item{\code{barcode}:}{Factor, the cell barcode for each molecule.
#' The levels contain the universe of all known barcodes.}
#' \item{\code{umi}:}{Integer, the processed UMI barcode in 2-bit encoding.} 
#' \item{\code{gem.group}:}{Integer, the GEM group.}
#' \item{\code{feature}:}{Integer, the index of the feature to which the molecule was assigned.
#' This refers to an entry in the \code{features} DataFrame, see below.}
#' \item{\code{count}:}{Integer, the number of reads mapped to this molecule.}
#' \item{\code{library}:}{Integer, the library index in cases where multiple libraries are present in the same file.
#' Only reported when \code{version="3"}.}
#' }
#' A field will not be present in the DataFrame if the corresponding \code{get.*} argument is \code{FALSE}, 
#' 
#' The second element of the list is \code{features}, a \link[S4Vectors]{DataFrame} where each row corresponds to a genomic feature.
#' This contains the following fields:
#' \describe{
#' \item{\code{ID}:}{Character, the identifier for the feature.
#' For genes, this is typically the Ensembl ID.}
#' \item{\code{Type}:}{Character, the type of the feature, e.g., \code{"Gene Expression"}, \code{"Antibody Capture"}.}
#' }
#' This is indexed by the \code{Feature} field in the \code{data} DataFrame.
#'
#' If \code{extract.library.info=TRUE}, the list will contain an additional element named \code{libraries}. 
#' This is a list of lists containing per-library information such as the \code{"library_type"}.
#' The \code{library} field in the \code{data} DataFrame indexes this list.
#' 
#' @details
#' Molecules that were not assigned to any gene have their \code{Feature} set to \code{nrow(features) + 1}.
#' These are removed when \code{keep.unmapped = FALSE}.
#' 
#' CellRanger 3.0 introduced a major change in the format of the molecule information files.
#' When \code{version="auto"}, the function will attempt to determine the version format of the file.
#' This can also be user-specified by setting \code{version} explicitly.
#' 
#' For files produced by version 2.2 of the CellRanger software, the length of the cell barcode is not given.
#' Instead, the barcode length is automatically inferred if \code{barcode.length=NULL} and \code{version="2"}.
#' Currently, version 1 of the 10X chemistry uses 14 nt barcodes, while version 2 uses 16 nt barcodes.
#' 
#' Setting any of the \code{get.*} arguments will (generally) avoid extraction of the corresponding field.
#' This can improve efficiency if that field is not necessary for further analysis.
#' Aside from the missing field, the results are guaranteed to be identical, i.e., same order and number of rows.
#' 
#' @author
#' Aaron Lun,
#' based on code by Jonathan Griffiths
#' 
#' @seealso
#' \code{\link{makeCountMatrix}}, which creates a count matrix from this information. 
#' 
#' @examples
#' # Mocking up some 10X HDF5-formatted data.
#' sim <- simulateMoleculeInformation()
#' temp <- tempfile(fileext=".h5")
#' writeMoleculeInformation(temp, sim)
#' 
#' # Reading the resulting file.
#' readMoleculeInformation(temp)
#' 
#' @references
#' Zheng GX, Terry JM, Belgrader P, and others (2017).
#' Massively parallel digital transcriptional profiling of single cells. 
#' \emph{Nat Commun} 8:14049.
#' 
#' 10X Genomics (2017).
#' Molecule info.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/2.2/output/molecule_info}
#' 
#' 10X Genomics (2018).
#' Molecule info.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info}
#'
#' @export
#' @importFrom rhdf5 h5read H5Fopen H5Fclose H5Dopen H5Dclose H5Dget_space H5Sget_simple_extent_dims H5Sclose
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
#' @importFrom jsonlite fromJSON
readMoleculeInformation <- function(
    path,
    barcode.length = NULL,
    keep.unmapped = FALSE, 
    get.barcode = TRUE,
    get.umi = TRUE,
    get.gem = TRUE,
    get.feature = TRUE,
    get.count = TRUE,
    get.library = TRUE,
    extract.library.info = TRUE,
    version = NULL
) {
    path <- path.expand(path) # protect against unexpanded tilde's.

    if (!is.null(version)) {
        version <- match.arg(version, c("2", "3"))
    } else {
        available <- h5ls(path, recursive=FALSE)
        if ("barcode_idx" %in% available$name) {
            version <- "3" 
        } else {
            version <- "2"
        }
    }

    data <- list()

    if (get.barcode) {
        if (version=="3") {
            all.barcodes <- as.vector(h5read(path, "barcodes"))
            barcodes.idx <- as.vector(h5read(path, "barcode_idx")) + 1L
            data$barcode <- factor(all.barcodes[barcodes.idx], all.barcodes)
        } else {
            data$barcode <- factor(read_cell_barcodes(path, "barcode", barcode.length))
        }
    }

    if (get.umi) {
        data$umi <- as.vector(h5read(path, "umi"))
    }

    if (get.gem) {
        data$gem.group <- as.vector(h5read(path, "gem_group"))
    }

    if (get.feature || !keep.unmapped) {
        # Both of these are zero-indexed by default.
        if (version == "3") {
            dname <- "feature_idx"
        } else {
            dname <- "gene"
        } 
        data$feature <- as.vector(h5read(path, dname)) + 1L 
    }

    if (get.count) {
        if (version=="3") {
            dname <- "count"
        } else {
            dname <- "reads"
        }
        data$count <- as.vector(h5read(path, dname))
    }

    if (version == "3" && get.library) {
        data$library <- as.vector(h5read(path, "library_idx")) + 1L
    }

    if (length(data) == 0) {
        # Just to ensure we get the right number of rows,
        # if there were no other fields requested.
        data <- (function() {
            fhandle <- H5Fopen(path)
            on.exit(H5Fclose(fhandle))
            dhandle <- H5Dopen(fhandle, "umi")
            on.exit(H5Dclose(dhandle), add=TRUE, after=FALSE)
            space <- H5Dget_space(dhandle)
            on.exit(H5Sclose(space), add=TRUE, after=FALSE)

            N <- H5Sget_simple_extent_dims(space)
            stopifnot(N$rank == 1L)
            make_zero_col_DFrame(N$size)
        })()
    } else {
        data <- do.call(DataFrame, data)
    }

    # Defining the set of all genes, removing unassigned gene entries.
    if (version=="3") {
        dname <- "features/id"
    } else {
        dname <- "gene_ids"
    }
    feat.ids <- as.vector(h5read(path, dname))

    if (!keep.unmapped) {
        keep <- (data$feature <= length(feat.ids))
        if (!get.feature) {
            data$feature <- NULL
        }
        data <- data[keep,,drop=FALSE]
    }

    # Don't define the total cell pool here, as higher level functions may want to use gem_group.
    output <- list(molecules = data, features = DataFrame(id = feat.ids))

    if (version == '3') {
        output$features$type <- as.vector(h5read(path, "features/feature_type"))

        if (extract.library.info) {
            lib.info <- as.vector(h5read(path, "library_info"))
            output$libraries <- fromJSON(lib.info, simplifyVector=FALSE)
        }
    }

    output
}
