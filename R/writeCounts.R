#' Write count data in the 10x format
#'
#' Create a directory containing the count matrix and cell/gene annotation from a sparse matrix of UMI counts, 
#' in the format produced by the CellRanger software suite.
#' 
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} containing a matrix of UMI counts.
#' @param path A string containing the path to the output directory (for \code{type="sparse"}) or file (for \code{type="HDF5"}).
#' @param assay Integer or string specifying the assay of \code{x} containing the count matrix to be written.
#' @param barcode.field String containing the name of the \code{\link[SummarizedExperiment]{colData}} column containing the cell barcodes.
#' If \code{NULL}, the column names of \code{x} are assumed to contain the barcodes.
#' @param feature.id.field String containing the name of the \code{\link[SummarizedExperiment]{rowData}} column containing the feature IDs.
#' If \code{NULL}, the row names of \code{x} are assumed to contain the feature IDs.
#' @param feature.name.field String containing the name of the \code{\link[SummarizedExperiment]{rowData}} column containing the feature names.
#' If \code{NULL}, the row names of \code{x} are assumed to contain the feature names.
#' @param feature.type.field String containing the name of the \code{\link[SummarizedExperiment]{rowData}} column containing the feature types.
#' If this name is not present in the \code{rowData}, the type defaults to \code{"Gene Expression"} for all rows of \code{x}.
#' Only used when \code{version="3"}.
#' @param overwrite A logical scalar specifying whether \code{path} should be overwritten if it already exists.
#' @param type String specifying the type of 10X format to save \code{x} to.
#' This is either a directory containing a sparse matrix with row/column annotation (\code{"sparse"})
#' or a HDF5 file containing the same information (\code{"HDF5"}).
#' @param genome String specifying the genome for storage when \code{type="HDF5"}.
#' This can be a character vector with one genome per feature if \code{version="3"}.
#' @param version String specifying the version of the CellRanger format to produce.
#' @param chemistry,original.gem.groups,library.ids 
#' Strings containing metadata attributes to be added to the HDF5 file for \code{type="HDF5"}.
#' Their interpretation is not formally documented and is left to the user's imagination.
#' @param num.threads Integer specifying the number of threads to use for counting the number of non-zero elements.
#' 
#' @details
#' This function will try to automatically detect the desired format based on whether \code{path} ends with \code{".h5"}.
#' If so, it assumes that \code{path} specifies a HDF5 file path and sets \code{type="HDF5"}.
#' Otherwise it will set \code{type="sparse"} under the assumption that \code{path} specifies a path to a directory.
#' 
#' Note that there were major changes in the output format for CellRanger version 3.0 to account for non-gene features such as antibody or CRISPR tags. 
#' Users can switch to this new format using \code{version="3"}.
#' See the documentation for \dQuote{latest} for this new format, otherwise see \dQuote{2.2} or earlier.
#'
#' The primary purpose of this function is to create files to use for testing \code{\link{read10xCounts}}.
#' In principle, it is possible to re-use the HDF5 matrices in \code{cellranger reanalyze}.
#' We recommend against doing so routinely due to CellRanger's dependence on undocumented metadata attributes that may change without notice.
#' 
#' @return 
#' For \code{type="sparse"}, a directory is produced at \code{path}.
#' If \code{version="2"}, this will contain the files \code{"matrix.mtx"}, \code{"barcodes.tsv"} and \code{"genes.tsv"}.
#' If \code{version="3"}, it will instead contain \code{"matrix.mtx.gz"}, \code{"barcodes.tsv.gz"} and \code{"features.tsv.gz"}.
#' 
#' For \code{type="HDF5"}, a HDF5 file is produced at \code{path} containing data in column-sparse format.
#' If \code{version="2"}, data are stored in the HDF5 group named \code{genome}.
#' If \code{version="3"}, data are stored in the group \code{"matrix"}.
#' 
#' A \code{TRUE} value is invisibly returned.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{read10xCounts}}, to read CellRanger matrices into R.
#' 
#' @examples
#' # Mocking up some count data.
#' library(Matrix)
#' my.counts <- matrix(rpois(1000, lambda=5), ncol=10, nrow=100)
#' my.counts <- as(my.counts, "dgCMatrix")
#' 
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(my.counts)
#' colnames(se) <- paste0("BARCODE-", seq_len(ncol(my.counts)))
#' rownames(se) <- paste0("ENSG0000", seq_len(nrow(my.counts)))
#' rowData(se)$Symbol <- paste0("GENE", seq_len(nrow(my.counts)))
#' rowData(se)$Type <- "Gene Expression"
#' 
#' # Writing this to file:
#' tmpdir <- tempfile()
#' writeCounts(tmpdir, se)
#' list.files(tmpdir)
#'
#' # Creating a version 3 HDF5 file:
#' tmph5 <- tempfile(fileext=".h5")
#' writeCounts(tmph5, se, version='3')
#' 
#' @references
#' 10X Genomics (2017).
#' Gene-Barcode Matrices.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/2.2/output/matrices}
#' 
#' 10X Genomics (2018).
#' Feature-Barcode Matrices.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices}
#' 
#' 10X Genomics (2018).
#' HDF5 Gene-Barcode Matrix Format.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/2.2/advanced/h5_matrices}
#' 
#' 10X Genomics (2018).
#' HDF5 Feature Barcode Matrix Format.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices}
#' 
#' @export
#' @importFrom SummarizedExperiment rowData colData
writeCounts <- function(
    path,
    x,
    assay = 1L,
    barcode.field = NULL, 
    feature.id.field = NULL,
    feature.name.field = "Symbol",
    feature.type.field = "Type",
    overwrite = FALSE,
    type = NULL,
    version = "3",
    genome = "unknown",
    compressed = (version == "3"),
    chemistry = "Single Cell 3' v3",
    original.gem.groups = 1L,
    library.ids = "custom",
    num.threads = 1
) {
    version <- match.arg(version, c("2", "3"))

    # Doing all the work on a temporary location next to 'path', as we have permissions there.
    # This avoids problems with 'path' already existing.
    temp.path <- tempfile(tmpdir=dirname(path)) 
    on.exit({ 
        if (file.exists(temp.path)) {
            unlink(temp.path, recursive=TRUE)
        } 
    })

    # Extracting all components.
    barcodes <- NULL
    if (is.null(barcode.field)) {
        barcodes <- colnames(x)
    } else {
        barcodes <- colData(x)[[barcode.field]]
    }
    if (is.null(barcodes)) {
        stop("no barcodes available at 'barcode.field'")
    }

    feature.ids <- NULL
    if (is.null(feature.id.field)) {
        feature.ids <- rownames(x)
    } else {
        feature.ids <- rowData(x)[[feature.id.field]]
    }
    if (is.null(feature.ids)) {
        stop("no feature IDs available at 'feature.id.field'")
    }

    feature.names <- NULL
    if (is.null(feature.name.field)) {
        feature.names <- rownames(x)
    } else {
        feature.names <- rowData(x)[[feature.name.field]]
    }
    if (is.null(feature.names)) {
        stop("no feature names available at 'feature.name.field'")
    }

    feature.types <- NULL
    if (version != "2") {
        if (feature.type.field %in% colnames(rowData(x))) {
            feature.types <- rowData(x)[[feature.type.field]]
        } else {
            feature.types <- rep("Gene Expression", nrow(x))
        }
    }

    x <- SummarizedExperiment::assay(x, assay)

    # Determining what format to save in.
    if (is.null(type)) {
        if (endsWith(path, ".h5")) {
            type <- "hdf5"
        } else {
            type <- "mtx"
        }
    }

    if (type == "mtx" || type == "prefix") {
        prefix <- NULL
        if (type == "prefix") {
            prefix <- basename(path)
        }

        .write_sparse(
            x = x,
            path = temp.path,
            barcodes = barcodes,
            feature.ids = feature.ids,
            feature.names = feature.names,
            feature.types = feature.types,
            version = version,
            compressed = compressed,
            prefix = prefix,
            num.threads = num.threads
        )
    } else {
        .write_hdf5(
            x = x,
            path = temp.path,
            genome = genome,
            barcodes = barcodes,
            feature.ids = feature.ids,
            feature.names = feature.names,
            feature.types = feature.types,
            version = version,
            chemistry = chemistry,
            original.gem.groups = original.gem.groups,
            library.ids = library.ids,
            num.threads = num.threads
        )
    }

    # We don't put this at the top as the write functions might fail; 
    # in which case, we would have deleted the existing 'path' for nothing.
    if (type == "prefix") {
        final.dir <- dirname(path)
        created <- list.files(temp.path)
        common <- intersect(list.files(final.dir), created)
        if (overwrite) {
            unlink(common, recursive=TRUE)
        } else if (length(common)) {
            stop("files with the specified 'path' prefix already exist")
        }
        file.rename(file.path(temp.path, created), file.path(final.dir, created))
    } else {
        if (overwrite) {
            unlink(path, recursive=TRUE)
        } else if (file.exists(path)) { 
            stop("specified 'path' already exists")
        }
        file.rename(temp.path, path)
    }

    invisible(NULL)
}

#' @importFrom utils write.table
#' @importFrom R.utils gzip
#' @importFrom DelayedArray type
#' @importFrom beachmat initializeCpp
.write_sparse <- function(
    path,
    x,
    barcodes,
    feature.ids,
    feature.names,
    feature.types,
    version,
    compressed,
    prefix,
    num.threads
) {
    gene.info <- data.frame(ID = feature.ids, Symbol = feature.names, stringsAsFactors=FALSE)

    bname <- "barcodes.tsv"
    mname <- "matrix.mtx"
    if (version == "3") {
        gene.info$Type <- feature.types
        fname <- "features.tsv"
    } else {
        fname <- "genes.tsv"
    }

    dir.create(path, recursive=TRUE, showWarnings=FALSE)
    if (!is.null(prefix)) {
        fname <- paste0(prefix, fname)
        bname <- paste0(prefix, bname)
        mname <- paste0(prefix, mname)
    }

    fpath <- file.path(path, fname)
    bpath <- file.path(path, bname)
    mpath <- file.path(path, mname)

    if (compressed) {
        mpath <- paste0(mpath, ".gz")
        bhandle <- gzfile(paste0(bpath, ".gz"), open="wb")
        fhandle <- gzfile(paste0(fpath, ".gz"), open="wb")
        on.exit({
            close(bhandle)
            close(fhandle)
        })
    } else {
        bhandle <- bpath
        fhandle <- fpath
    }

    write(barcodes, file=bhandle)
    write.table(gene.info, file=fhandle, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

    write_mm(
        initializeCpp(x),
        mpath,
        compressed = compressed,
        is_integer = (type(x) == "integer"),
        num_threads = num.threads
    )
}

#' @importFrom beachmat initializeCpp
#' @importFrom rhdf5 h5createFile h5createGroup h5write h5writeAttribute H5Gopen H5Fopen H5Gclose H5Fclose
#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
.write_hdf5 <- function(
    path,
    genome,
    x,
    barcodes,
    feature.ids,
    feature.names,
    feature.types,
    version,
    chemistry,
    original.gem.groups,
    library.ids,
    num.threads
) {
    path <- path.expand(path) # protect against tilde's.
    h5createFile(path)

    if (version=="3") {
        group <- "matrix"
    } else {
        group <- genome
    }
    h5createGroup(path, group)

    h5write(barcodes, file=path, name=paste0(group, "/barcodes"))

    # Saving feature information.
    if (version=="3") {
        h5createGroup(path, file.path(group, "features"))

        h5write(feature.ids, file=path, name=paste0(group, "/features/id"))
        h5write(feature.names, file=path, name=paste0(group, "/features/name"))
        h5write(rep(feature.types, length.out=length(feature.ids)),
            file=path, name=paste0(group, "/features/feature_type"))

        h5write("genome", file=path, name=paste0(group, "/features/_all_tag_keys"))
        h5write(rep(genome, length.out=length(feature.ids)),
            file=path, name=paste0(group, "/features/genome"))

        # Writing attributes.
        (function() {
            h5f <- H5Fopen(path)
            on.exit(H5Fclose(h5f), add=TRUE, after=FALSE)

            h5g <- H5Gopen(h5f, "/")
            on.exit(H5Gclose(h5g), add=TRUE, after=FALSE)

            h5writeAttribute(chemistry, h5obj=h5g, name="chemistry_description")
            h5writeAttribute("matrix", h5obj=h5g, name="filetype")
            h5writeAttribute(library.ids, h5obj=h5g, name="library_ids")
            h5writeAttribute(original.gem.groups, h5obj=h5g, name="original_gem_groups")
            h5writeAttribute(as.integer(version) - 1L, h5obj=h5g, name="version") # this is probably correct.
        })()

    } else {
        h5write(feature.ids, file=path, name=paste0(group, "/genes"))
        h5write(feature.names, file=path, name=paste0(group, "/gene_names"))
    }

    write_hdf5_counts(
        initializeCpp(x),
        path = path,
        group = group,
        num_threads = num.threads
    )
}
