#include "Rcpp.h"

#include <stdexcept>
#include <string>

#include "sanisizer/sanisizer.hpp"

#include "utils.h"

//[[Rcpp::export(rng=false)]]
Rcpp::IntegerVector encode_sequences(Rcpp::StringVector sequences) {
    const auto nseq = sequences.size();
    auto output = sanisizer::create<Rcpp::IntegerVector>(nseq);

    for (I<decltype(nseq)> i = 0; i < nseq; ++i) {
        Rcpp::String current(sequences[i]);
        const char* ptr = current.get_cstring();
        const auto len = Rf_length(current.get_sexp());
        if (len > 15) {
            throw std::runtime_error("32-bit signed integers cannot encode sequences longer than 15 nt");
        }

        int encoded = 0;
        for (I<decltype(len)> j = 0; j < len; ++j) {
            encoded <<= 2;
            switch (ptr[j]) {
                case 'A':
                    break;
                case 'C':
                    encoded += 1;
                    break;
                case 'G':
                    encoded += 2;
                    break;
                case 'T':
                    encoded += 3;
                    break;
                default:
                    throw std::runtime_error("unrecognized character in '" + std::string(ptr) + "'");
            }
        }

        output[i] = encoded;
    }

    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::StringVector decode_sequences(Rcpp::IntegerVector encoded, Rcpp::IntegerVector lengths) {
    const auto nseq = encoded.size();
    if (nseq != lengths.size()) {
        throw std::runtime_error("lengths of 'sequences' and 'lengths' should be the same");
    }
    auto output = sanisizer::create<Rcpp::StringVector>(nseq);

    std::string payload;
    const char* bases = "ACGT";
    for (I<decltype(nseq)> i = 0; i < nseq; ++i) {
        payload.clear();

        auto val = encoded[i];
        const auto len = lengths[i];
        for (I<decltype(len)> j = 0; j < len; ++j) {
            payload += bases[val & 0x3];
            val >>= 2;
        }

        std::reverse(payload.begin(), payload.end());
        output[i] = payload;
    }

    return output;
}
