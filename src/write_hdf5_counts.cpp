#include "Rtatami.h"
#include "Rcpp.h"

#include "H5Cpp.h"
#include "tatami/tatami.hpp"
#include "tatami_hdf5/tatami_hdf5.hpp"

#include <string>
#include <memory>
#include <cstdint>

#include "utils.h"

//[[Rcpp::export(rng=false)]]
SEXP write_hdf5_counts(Rcpp::RObject ptr, std::string path, std::string group, int num_threads) {
    Rtatami::BoundNumericPointer mptr(ptr);

    H5::H5File fhandle(path, H5F_ACC_RDWR);
    auto ghandle = fhandle.openGroup(group);

    tatami_hdf5::WriteCompressedSparseMatrixOptions opt;
    opt.columnar = tatami_hdf5::WriteStorageLayout::COLUMN;
    opt.data_type = tatami_hdf5::WriteStorageType::UINT32;
    opt.index_type = tatami_hdf5::WriteStorageType::UINT32;
    // opt.indptr_type = tatami_hdf5::WriteStorageType::UINT32;
    opt.num_threads = num_threads;

    const auto& mat = *(mptr->ptr);
    tatami_hdf5::write_compressed_sparse_matrix(mat, ghandle, opt);

    // Adding the shape, while we're here.
    {
        hsize_t ndim = 2;
        H5::DataSpace dspace(1, &ndim);
        auto dhandle = ghandle.createDataSet("shape", H5::PredType::NATIVE_UINT64, dspace);

        std::array<std::uint64_t, 2> dims { 
            sanisizer::cast<std::uint64_t>(mat.nrow()),
            sanisizer::cast<std::uint64_t>(mat.ncol())
        };
        dhandle.write(dims.data(), H5::PredType::NATIVE_UINT64, dspace);
    }

    return R_NilValue;
}
