#' Load data from a 10X Genomics experiment
#' 
#' Creates a \link[SingleCellExperiment]{SingleCellExperiment} from the CellRanger output directories for 10X Genomics data.
#' 
#' @param samples A list of sample information.
#' Each element corresponds to a sample and is itself a list returned by \code{configureSampleForReadCounts}.
#'
#' Alternatively, a character vector of paths.
#' Each entry will be interpreted as described below for \code{path} with an automatically-detected \code{type}.
#'
#' In both cases, the list or vector can be named, in which case the names are used to identify each sample in the output object.
#' If unnamed, the sample identity is set to the path instead.
#' @param column.names Boolean indicating whether the columns of the output object should be named with the cell barcodes.
#' @param row.names String specifying whether to use Ensembl IDs ("ID") or gene symbols ("Symbol") as row names. 
#' For symbols, the Ensembl ID will be appended to disambiguate rows where the same symbol corresponds to multiple Ensembl IDs.
#' @param delayed Boolean indicating whether sparse matrices should be wrapped in \link[DelayedArray]{DelayedArray}s before combining.
#' This saves memory by avoiding the creation of a new matrix; it also avoids integer overflow for \link[Matrix]{dgCMatrix} instances.
#' Only relevant for multiple \code{samples}.
#' @param intersect.rows Boolean indicating whether to take the intersection of common row names across all samples.
#' If \code{FALSE}, differences in row names across samples will cause an error to be raised.
#' @param BPPARAM A \link[BiocParallel]{BiocParallelParam} object specifying how loading should be parallelized for multiple \code{samples}.
#' @param path String specifying the path to the CellRanger output for a single sample.
#' \itemize{
#' \item If \code{type="hdf5"}, \code{path} should be a HDF5 file in the 10X sparse HDF5 format. 
#' Users may need to set \code{genome} if it cannot be automatically determined when \code{version="2"}.
#' \item If \code{type="mtx"}, \code{path} should refer to a directory containing:
#' \itemize{
#' \item \code{"matrix.mtx"}, a Matrix Market file containing a feature-by-barcode UMI count matrix.
#' \item \code{"barcodes.tsv"}, a tab-delimited file containing per-barcode information.
#' \item \code{"features.tsv"}, a tab-delimited file containing information for each genomic feature (usually genes).
#' For \code{version="2"}, this will be called \code{"genes.tsv"} instead.
#' }
#' Any of these files may be Gzip-compressed, in which case they should have an additional \code{".gz"} file extension.
#' \item If \code{type="prefix"}, \code{path} should be a path prefix for all files associated with this sample.
#' The files should follow the same format described above for \code{type="mtx"}.
#' For example, if \code{path="xyz_"}, the files associated with this sample should be \code{"xyz_matrix.mtx"}, \code{"xyz_barcodes.tsv"}, etc.
#' Again, any of these files may be Gzip-compressed, in which case they should be suffixed with \code{".gz"}.
#' }
#' @param type String specifying the type of 10X format to read data from.
#' If \code{NULL}, the function will attempt to automatically detect this from the specified files -
#' \code{"hdf5"} if the path ends with \code{".h5"}, \code{"mtx"} if the path refers to a directory, and \code{"prefix"} otherwise.
#' @param version String specifying the version of the 10X format to read data from.
#' If \code{NULL}, the function will attempt to automatically detect this from the specified files.
#' For \code{type="mtx"}, this is based on whether the per-feature information is stored in a \code{"features.tsv*"} or \code{"genes.tsv*"} file.
#' For \code{type="hdf5"}, this is based on whether there is a top-level \code{"matrix"} group with a \code{"matrix/features"} subgroup in the file.
#' @param genome String specifying the genome, i.e., the HDF5 group containing the data.
#' If \code{NULL}, the function will attempt to automatically detect this from the specified files.
#' Only used if \code{type="hdf5"} and \code{version='2'}.
#' @param compressed Boolean indicating whether the text files are compressed for \code{type="mtx"} or \code{"prefix"}.
#' If \code{NULL}, the function will attempt to automatically detect this from the specified files, i.e., whether the file paths end with \code{".gz"}. 
#' @param mtx.two.pass Boolean indicating whether to use a two-pass approach for loading data from a Matrix Market file.
#' This reduces peak memory usage at the cost of some additional runtime. 
#' Only relevant when \code{type="mtx"} or \code{type="prefix"}.
#' @param mtx.class String specifying the class of the output matrix when \code{type="mtx"} or \code{type="prefix"}.
#' @param mtx.threads Boolean specifying the number of threads to use for reading Matrix Market files.
#' Only relevant when \code{type="mtx"} or \code{type="prefix"}.
#' 
#' @return 
#' For \code{readCounts}, a \link[SingleCellExperiment]{SingleCellExperiment} object containing count data for each genomic feature (row) and cell (column) across all \code{samples}.
#' \itemize{
#' \item The \code{\link[SummarizedExperiment]{rowData}} contains the \code{"ID"} and \code{"Symbol"} columns.
#' The former is the gene identifier (usually Ensembl), while the latter is the gene name.
#' If \code{version="3"}, it will also contain the \code{"Type"} field specifying the type of feature (e.g., gene, antibody tag, CRISPR guide).
#' \item The \code{\link[SummarizedExperiment]{rowRanges}} contains the genomic coordinates for each gene.
#' Only provided for \code{version="3"}, 
#' \item The \code{\link[SummarizedExperiment]{colData}} will contain the \code{"Sample"} and \code{"Barcode"} columns.
#' The former contains the name of the sample (or if not supplied, the path in \code{samples}) from which each column was obtained.
#' The latter contains the cell barcode sequence and GEM group for each cell library. 
#' \item Rows are named with the gene identifier.
#' If \code{column.names=TRUE}, each column is named by the cell barcode if there is only a single sample.
#' For multiple samples, the index of each sample in \code{samples} is concatenated to the cell barcode to form the column name.
#' This avoids problems with multiple instances of the same cell barcodes in different samples.
#' \item The \code{\link[SummarizedExperiment]{assays}} will contain a single \code{"counts"} matrix, containing UMI counts for each gene in each cell in each sample.
#' Matrices are combined by column if multiple \code{samples} were specified.
#' Note that the matrix representation will depend on the format of the \code{samples}:
#' \itemize{
#' \item For \code{type="prefix"} or \code{type="mtx"}, the class of the matrix is determined by \code{mtx.class}.
#' Otherwise, for \code{type="hdf5"}, the matrix is a \link[HDF5Array]{HDF5Matrix}.
#' \item If multiple samples are present, the class of the matrix is determined by the return type of \code{cbind} on the per-sample matrices.
#' This can be forced to be a \link[DelayedArray]{DelayedMatrix} by setting \code{delayed=TRUE}. 
#' }
#' \item The \code{\link[S4Vectors]{metadata}} contains a \code{"Samples"} field, containing the input \code{samples} character vector.
#' }
#'
#' For \code{configureSampleForReadCounts}, a list of configuration parameters for each sample.
#' This should be passed to \code{samples} in \code{readCounts}.
#' 
#' @details
#' This function has a long and storied past.
#' It was originally developed as the \code{read10xResults} function in \pkg{scater}, inspired by the \code{Read10X} function from the \pkg{Seurat} package.
#' It was then migrated to \pkg{DropletUtils} in an effort to consolidate some 10X-related functionality across various packages.
#' Finally, it was moved to \pkg{crio} to separate the input/output utilities from analysis functions.
#'
#' Note that user-level manipulation of sparse matrices requires loading of the \pkg{Matrix} package.
#' Otherwise, calculation of \code{rowSums}, \code{colSums}, etc. will result in errors.
#' 
#' @author
#' Davis McCarthy, with modifications from Aaron Lun.
#' 
#' @seealso
#' \code{\link[SingleCellExperiment]{splitAltExps}}, to split alternative feature sets (e.g., antibody tags) into their own Experiments.
#' 
#' \code{\link{write10xCounts}}, to create 10X-formatted file(s) from a count matrix.
#' 
#' @examples
#' # Mocking up some 10X genomics output.
#' example(write10xCounts, echo=FALSE)
#' 
#' # Reading it in.
#' sce10x <- readCounts(tmpdir)
#' 
#' # Works with multiple samples.
#' sce10x2 <- readCounts(c(tmpdir, tmpdir))
#' 
#' @references
#' Zheng GX, Terry JM, Belgrader P, and others (2017).
#' Massively parallel digital transcriptional profiling of single cells. 
#' \emph{Nat Commun} 8:14049.
#' 
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
#' @importFrom S4Vectors cbind DataFrame ROWNAMES<- ROWNAMES extractROWS
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom GenomicRanges GRanges
#' @importFrom DelayedArray DelayedArray
readCounts <- function(
    samples, 
    column.names = FALSE, 
    row.names = c("id", "symbol"),
    delayed = FALSE,
    intersect.rows = FALSE,
    BPPARAM= SerialParam()
) {
    row.names <- match.arg(row.names)

    if (is.character(samples)) {
        old.names <- names(samples)
        samples <- lapply(samples, configureSampleForReadCounts)
        names(samples) <- old.names
    } else if (inherits(samples, "readCountsSample")) {
        samples <- list(samples)
    }

    sample.names <- names(samples)
    sample.paths <- vapply(samples, function(x) x$path, FUN.VALUE="")
    if (is.null(sample.names)) {
        sample.names <- sample.paths
    }

    load.out <- bplapply(samples, FUN=.tenx_loader, BPPARAM=BPPARAM)

    nsets <- length(samples)
    full_data <- vector("list", nsets)
    gene_info_list <- vector("list", nsets)
    cell_info_list <- vector("list", nsets)

    for (i in seq_len(nsets)) { 
        current <- load.out[[i]]
        full_data[[i]] <- current$mat

        rr <- current$gene.info
        if (row.names == "id") {
            rns <- rr$ID
        } else {
            rns <- rr$Symbol
            if (anyDuplicated(rns)) {
                dup.name <- rns %in% rns[duplicated(rns)]
                rns[dup.name] <- paste0(rns[dup.name], "_", rr$ID[dup.name])
            }
        }

        ROWNAMES(rr) <- rns
        gene_info_list[[i]] <- rr

        cell.names <- current$cell.names
        cell_info_list[[i]] <- DataFrame(
            Sample = rep(sample.names[i], length(cell.names)), 
            Barcode = cell.names,
            row.names=NULL
        )
    }

    # Checking gene correctness and handling differences, or failing.
    gene_info <- gene_info_list[[1]]
    if (!intersect.rows) {
        if (nsets > 1 && length(unique(gene_info_list)) != 1L) {
            stop("gene information differs between runs")
        }
    } else {
        all.genes <- lapply(gene_info_list, ROWNAMES)
        if (length(unique(all.genes)) > 1) {
            common.genes <- Reduce(intersect, all.genes)
            for (i in seq_along(full_data)) {
                m <- match(common.genes, all.genes[[i]])
                full_data[[i]] <- full_data[[i]][m,,drop=FALSE]
            }
            gene_info <- extractROWS(gene_info, common.genes)
        }
    }

    # Forming the full data matrix.
    if (length(full_data) > 1) {
        if (delayed) {
            full_data <- lapply(full_data, DelayedArray)
        }
        full_data <- do.call(cbind, full_data)
    } else {
        full_data <- full_data[[1]]
    }

    # Adding the cell data.
    cell_info <- do.call(rbind, cell_info_list)
    if (column.names) {
        if (nsets == 1L) {
            cnames <- cell_info$Barcode
        } else {
            sid <- rep(seq_along(cell_info_list), vapply(cell_info_list, nrow, 1L))
            cnames <- paste0(sid, "_", cell_info$Barcode)
        }
        colnames(full_data) <- cnames
    }

    SingleCellExperiment(
        list(counts = full_data),
        rowData = gene_info,
        colData = cell_info,
        metadata=list(Samples=sample.paths)
    )
}

.type_chooser2 <- function(path, type) {
    if (!is.null(type)) {
        match.arg(type, c("mtx", "hdf5", "prefix"))
    } else if (grepl("\\.h5", path)) {
        "hdf5"
    } else if (dir.exists(path)) {
        "mtx"
    } else {
        "prefix"
    }
}

#' @export
#' @rdname readCounts
configureSampleForReadCounts <- function(
    path,
    type = NULL,
    version = NULL,
    genome = NULL, 
    compressed = NULL, 
    mtx.two.pass = FALSE,
    mtx.class = c("CsparseMatrix", "SVT_SparseMatrix"),
    mtx.threads = 1
) {
    output <- list(
        path = path,
        type = .type_chooser2(path, type),
        version = version,
        genome = genome,
        compressed = compressed,
        mtx.two.pass = mtx.two.pass,
        mtx.class = match.arg(mtx.class),
        mtx.threads = mtx.threads
    )
    class(output) <- "readCountsSample"
    output
}

.tenx_loader <- function(sample) {
    cur.type <- sample$type
    if (cur.type=="mtx" || cur.type == "prefix") {
        .read_from_sparse(
            sample$path,
            version = sample$version,
            is.prefix = (cur.type == "prefix"),
            compressed = sample$compressed,
            mtx.two.pass = sample$mtx.two.pass,
            mtx.class = sample$mtx.class,
            mtx.threads = sample$mtx.threads
        )
    } else {
        .read_from_hdf5(
            sample$path,
            genome = sample$genome,
            version = sample$version
        )
    }
}

#' @importFrom methods as new
#' @importFrom utils read.delim head
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols<- DataFrame
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom SparseArray SVT_SparseMatrix 
.read_from_sparse <- function(path, version, is.prefix, compressed, mtx.two.pass, mtx.class, mtx.threads) {
    if (is.prefix) {
        FUN <- paste0
    } else {
        FUN <- file.path
    }

    if (is.null(version)) {
        target <- FUN(path, "features.tsv")
        target <- .check_for_compressed(target, compressed, error=FALSE)
        if (file.exists(target)) {
            version <- "3"
        } else {
            version <- "2"
        }
    }

    bname <- "barcodes.tsv"
    mname <- "matrix.mtx"
    if (version=="3") {
        gname <- "features.tsv"
    } else {
        gname <- "genes.tsv"
    }

    barcode.loc <- FUN(path, bname)
    gene.loc <- FUN(path, gname)
    matrix.loc <- FUN(path, mname)

    barcode.loc <- .check_for_compressed(barcode.loc, compressed)
    gene.loc <- .check_for_compressed(gene.loc, compressed)
    matrix.loc <- .check_for_compressed(matrix.loc, compressed)

    gene.info <- read.delim(gene.loc, header=FALSE, stringsAsFactors=FALSE, quote="", comment.char="")
    possible.names <- c("ID", "Symbol", "Type", "Chromosome", "Start", "End")
    colnames(gene.info) <- head(possible.names, ncol(gene.info))

    if (ncol(gene.info) > 3) {
        # Default ARC-seq reference seems to give empty names for mitochondrial genes.
        # Hey, I don't make the rules.
        keep <- gene.info$Chromosome == "" & grepl("^MT-", gene.info$Symbol)
        if (any(keep)) {
            gene.info[keep,"Chromosome"] <- "chrM"
        }

        # Newer cellranger versions like to add coordinates, so 
        # let's throw it into the GRanges for fun.
        gr <- GRanges(gene.info$Chromosome, IRanges(gene.info$Start, gene.info$End))
        mcols(gr) <- DataFrame(gene.info[,1:3])
        gene.info <- gr 
    }

    raw_mat <- read_mm(matrix.loc, two_pass=mtx.two.pass, class_name=mtx.class, threads=mtx.threads)
    if (mtx.class == "CsparseMatrix") {
        # Don't use sparseMatrix as this seems to do an unnecessary roundtrip through the triplet form.
        mat <- new("dgCMatrix", Dim=raw_mat$dim, i=raw_mat$contents$i, x=raw_mat$contents$x, p=raw_mat$contents$p) 
    } else {
        mat <- new("SVT_SparseMatrix", SVT=raw_mat$contents$list, dim=raw_mat$dim, type=raw_mat$contents$type)
    }

    list(
        mat=mat,
        cell.names=readLines(barcode.loc),
        gene.info=gene.info
    )
}

.check_for_compressed <- function(path, compressed, error=TRUE) {
    original <- path
    if (isTRUE(compressed)) {
        path <- paste0(path, ".gz")
    } else if (is.null(compressed) && !file.exists(path)) {
        path <- paste0(path, ".gz")
        if (error && !file.exists(path)) {
            # Explicit error here to avoid users getting confused.
            stop(sprintf("cannot find '%s' or its gzip-compressed form", original)) 
        }
    }
    path
}

#' @importFrom rhdf5 h5ls h5read
#' @importFrom HDF5Array TENxMatrix
#' @importFrom utils head
.read_from_hdf5 <- function(path, genome=NULL, version) {
    path <- path.expand(path) # eliminate tilde's.
    available <- h5ls(path, recursive=FALSE)
    available <- available[available$otype=="H5I_GROUP",]

    if (is.null(version)) {
        if ("matrix" %in% available$name) {
            version <- "3"
        } else {
            version <- "2"
        }
    }

    if (version=="2") {
        group <- genome
        if (is.null(group)) {
            if (nrow(available) > 1L) {
                to.see <- head(available$name, 3)
                if (length(to.see)==3L) {
                    to.see[3] <- "..."
                }
                stop("more than one available group (", paste(to.see, collapse=", "), ")")
            } else if (nrow(available) == 0L) {
                stop("no available groups")
            }
            group <- available$name
        }

        gene.info <- data.frame(
            ID=as.character(h5read(path, paste0(group, "/genes"))),
            Symbol=as.character(h5read(path, paste0(group, "/gene_names"))),
            stringsAsFactors=FALSE
        )
    } else {
        group <- "matrix"
        gene.info <- data.frame(
            ID=as.character(h5read(path, paste0(group, "/features/id"))),
            Symbol=as.character(h5read(path, paste0(group, "/features/name"))),
            Type=as.character(h5read(path, paste0(group, "/features/feature_type"))),
            stringsAsFactors=FALSE
        )
    }

    mat <- TENxMatrix(path, group)
    dimnames(mat) <- NULL # for consistency.
    list(
        mat=mat,
        cell.names=as.character(h5read(path, paste0(group, "/barcodes"))),
        gene.info=gene.info
    )
}
