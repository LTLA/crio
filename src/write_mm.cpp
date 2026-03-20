#include "Rtatami.h"
#include "Rcpp.h"

#include "tatami/tatami.hpp"
#include "tatami_mtx/tatami_mtx.hpp"
#include "byteme/byteme.hpp"

#include <string>
#include <memory>

#include "utils.h"

//[[Rcpp::export(rng=false)]]
SEXP write_mm(Rcpp::RObject ptr, std::string path, bool compressed, bool is_integer, int num_threads) {
    Rtatami::BoundNumericPointer mptr(ptr);

    std::unique_ptr<byteme::Writer> writer;
    if (compressed) {
        writer.reset(new byteme::GzipFileWriter(path.c_str(), {}));
    } else {
        writer.reset(new byteme::RawFileWriter(path.c_str(), {}));
    }

    // Writing our own banner.
    writer->write("%%MatrixMarket matrix coordinate");
    if (is_integer) {
        writer->write(" integer");
    } else {
        writer->write(" real");
    }
    writer->write(" general\n");

    tatami_mtx::WriteMatrixOptions opt;
    opt.coordinate = true;
    opt.banner = false;
    opt.num_threads = num_threads;
    if (is_integer) {
        opt.format = std::chars_format::fixed;
    }

    tatami_mtx::write_matrix(*(mptr->ptr), *writer, opt);
    writer->finish();
    return R_NilValue;
}
