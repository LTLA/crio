#' Make a count matrix from molecule-level information
#'
#' Construct a count matrix from per-molecule information, typically the barcode and feature to which each molecule was assigned.
#'
#' @param num.features Integer specifying the total number of features used in the reference annotation.
#' @param feature Integer vector of length equal to the number of molecules, containing the feature to which molecule was assigned.
#' A value greater than \code{num.features} indicates that the corresponding molecule is unassigned.
#' @param barcode Factor of length equal to the number of molecules, containing the cell barcode to which each molecule was assigned.
#' Levels are usualy barcode sequences.
#' @param counts Integer vector of length equal to the number of molecules, containing the read counts for each molecule.
#' If supplied, the count matrix refers to the number of reads assigned to each gene in each barcode, not the number of molecules (i.e., UMIs).
#' @param num.threads Integer specifying the number of threads to use.
#' @param class String specifying the output class.
#'
#' @return A sparse matrix of the specified \code{class}.
#' Each row is a feature and each column is a barcode.
#' Each entry is the number of molecules (if \code{counts=NULL}) or reads (otherwise) for the corresponding feature and barcode.
#' Columns are named after the corresponding levels of \code{barcode}.
#'
#' @author Aaron Lun
#'
#' @examples
#' # Mocking up some 10X molecule information.
#' sim <- simulateMoleculeInformation(num.features=20)
#'
#' # Formatting it into a matrix.
#' countMolecules(20, sim$molecules$feature, sim$molecules$barcode)
#'
#' @seealso
#' \code{\link{readMolecules}}, which is typically used to load molecule information from the HDF5 file produced by CellRanger.
#'
#' @export
countMolecules <- function(
    num.features,
    feature,
    barcode,
    count = NULL,
    class = c("dgCMatrix", "SVT_SparseMatrix"),
    num.threads = 1
) {
    barcode <- as.factor(barcode)
    use.svt <- match.arg(class) == "SVT_SparseMatrix"
    out <- make_count_matrix(
        num.features,
        nlevels(barcode),
        feature,
        as.integer(barcode),
        count,
        use.svt,
        num.threads
    )

    if (use.svt) {
        new("SVT_SparseMatrix",
            SVT=out,
            dim=c(as.integer(num.features), nlevels(barcode)),
            dimnames=list(NULL, levels(barcode)),
            type="integer"
        )

    } else {
        new("dgCMatrix", 
            i=out[[1]],
            x=out[[2]],
            p=out[[3]],
            Dim=c(as.integer(num.features), nlevels(barcode)),
            Dimnames=list(NULL, levels(barcode))
        )
    }
}
