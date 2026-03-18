#include "Rcpp.h"

#include "H5Cpp.h"

#include <stdexcept>
#include <cstdint>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

//[[Rcpp::export(rng=false)]]
Rcpp::StringVector read_cell_barcodes(std::string fname, std::string dname, Rcpp::Nullable<Rcpp::IntegerVector> barcodelen) {
    H5::H5File h5file(fname.c_str(), H5F_ACC_RDONLY);
    H5::DataSet h5data = h5file.openDataSet(dname.c_str());

    if (h5data.getTypeClass() != H5T_INTEGER) {
        throw std::runtime_error("cell barcodes should be encoded as integers");
    }

    H5::DataSpace dataspace = h5data.getSpace();
    if (dataspace.getSimpleExtentNdims()!=1) {
        throw std::runtime_error("cell barcodes should be a one-dimensional array");
    }
    hsize_t dims_out;
    dataspace.getSimpleExtentDims(&dims_out, NULL);

    H5::DataSpace memspace(1, &dims_out);
    memspace.selectAll();
    dataspace.selectAll();
    std::vector<std::uint64_t> encoded(dims_out);
    h5data.read(encoded.data(), H5::PredType::NATIVE_UINT64, memspace, dataspace);
   
    // Guessing the barcode length. 
    int blen = 0;
    if (barcodelen.isNull()) {
        if (encoded.size()) {
            blen = std::ceil(std::log2(*std::max_element(encoded.begin(), encoded.end())) * 0.5); // divide by 2 as we really want log4.
        }
    } else {
        Rcpp::IntegerVector blens(barcodelen); 
        if (blens.size() != 1) {
            throw std::runtime_error("barcode length should be an integer scalar");
        }
        blen = blens[0];
        if (blen < 0) {
            throw std::runtime_error("barcode length should be a non-negative integer");
        }
    }

    Rcpp::StringVector output(dims_out);
    auto oIt = output.begin();
    std::vector<char> ref(blen + 1, '\0');   
    const char* bases = "ACGT";

    for (auto enc : encoded) {
        for (int pos = 0; pos < blen; ++pos) {
            ref[blen - pos - 1] = bases[enc & 0x3];
            enc >>= 2;
        }

        (*oIt) = Rcpp::String(ref.data());
        ++oIt;
    }

    return output;    
}

//[[Rcpp::export(rng=false)]]
SEXP write_cell_barcodes(Rcpp::StringVector barcodes, std::string fname, std::string dname) {
    const auto num_barcodes = barcodes.size();
    std::vector<std::uint64_t> encoded;
    encoded.reserve(num_barcodes); 

    int len = -1;
    for (auto bc0 : barcodes) {
        Rcpp::String bc(bc0);
        const char* ptr = bc.get_cstring();
        const int curlen = Rf_length(bc.get_sexp());
        if (len < 0) {
            len = curlen;
            if (len > 32) {
                throw std::runtime_error("barcode is too long to encode in a 64-bit integer");
            }
        } else if (len != curlen) {
            throw std::runtime_error("all barcodes should have the same length");
        }

        std::uint64_t out = 0;
        for (int i = 0; i < len; ++i) {
            out <<= 2;
            switch (ptr[i]) {
                case 'A': case 'a':
                    break;
                case 'C': case 'c':
                    out += 1;
                    break;
                case 'G': case 'g':
                    out += 2;
                    break;
                case 'T': case 't':
                    out += 3;
                    break;
                default:
                    throw std::runtime_error("unrecognized base '" + std::string(1, ptr[i]) + "' in the barcode");
            }
        }

        encoded.push_back(out);
    }

    hsize_t length = num_barcodes;
    H5::DataSpace dspace(1, &length);
 	H5::DSetCreatPropList plist;

    if (length) {
        plist.setDeflate(6);
        hsize_t chunk = std::max(static_cast<hsize_t>(1000), static_cast<hsize_t>(std::sqrt(num_barcodes)));
        if (chunk > length) {
            plist.setChunk(1, &length);
        } else {
            plist.setChunk(1, &chunk);
        }
    }

    H5::H5File fhandle(fname.c_str(), H5F_ACC_RDWR);
    auto dhandle = fhandle.createDataSet(dname, H5::PredType::NATIVE_UINT64, dspace, plist);
    dhandle.write(encoded.data(), H5::PredType::NATIVE_UINT64, dspace);

    return R_NilValue;
}
