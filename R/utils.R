#' @import rhdf5
h5write_type <- function(x, path, group, name, type) {
    fhandle <- H5Fopen(path, "H5F_ACC_RDWR")
    on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)

    if (!is.null(group)) {
        ghandle <- H5Gopen(fhandle, group)
        on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)
    } else {
        ghandle <- fhandle
    }

    dcpl <- H5Pcreate("H5P_DATASET_CREATE")
    on.exit(H5Pclose(dcpl), add = TRUE, after=FALSE)

    H5Pset_chunk(dcpl, max(1000, sqrt(length(x))))
    H5Pset_shuffle(dcpl)
    H5Pset_deflate(dcpl, level = 6)

    space <- H5Screate_simple(length(x))
    on.exit(H5Sclose(space), add = TRUE, after=FALSE)

    dhandle <- H5Dcreate(ghandle, name, type, space, dcpl=dcpl)
    on.exit(H5Dclose(dhandle), add = TRUE, after = FALSE)

    H5Dwrite(dhandle, x)
}
