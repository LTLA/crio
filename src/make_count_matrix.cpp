#include "Rcpp.h"

#include "sanisizer/sanisizer.hpp"
#include "subpar/subpar.hpp"

#include <vector>
#include <stdexcept>
#include <type_traits>
#include <algorithm>

#include "utils.h"

template<typename Counts_>
Rcpp::RObject make_count_matrix_raw(int nrow, int ncol, Rcpp::IntegerVector features, Rcpp::IntegerVector barcodes, Counts_ counts, bool create_svt, int num_threads) { 
    const auto num_barcodes = barcodes.size();
    if (num_barcodes != features.size()) {
        throw std::runtime_error("'barcodes' and 'features' should be of the same length");
    }

    constexpr bool use_counts = !std::is_same<Counts_, bool>::value;
    if constexpr(use_counts) {
        if (num_barcodes != counts.size()) {
            throw std::runtime_error("'barcodes' and 'counts' should be of the same length");
        }
    }

    // Counting the number of elements for preallocation.
    typedef I<decltype(std::declval<Rcpp::IntegerVector>().size())> RLen;
    std::vector<RLen> cumulative_elements(sanisizer::sum<typename std::vector<RLen>::size_type>(ncol, 1));
    for (I<decltype(num_barcodes)> i = 0; i < num_barcodes; ++i) {
        // No need to check barcodes[i] <= ncol, as integer codes will be no greater than nlevels.
        if (barcodes[i] > 0 && features[i] > 0 && features[i] <= nrow) {
            // Barcode indices are 1-based but we have to add 1 anyway to get the right cumulative entry.
            ++cumulative_elements[barcodes[i]]; 
        }
    }

    for (int c = 1; c <= ncol; ++c) {
        cumulative_elements[c] += cumulative_elements[c - 1];
    }

    // Populating the vector of per-column features.
    auto offsets = cumulative_elements;
    auto re_features = [&](){
        if constexpr(use_counts) {
            return sanisizer::create<std::vector<std::pair<int, int> > >(cumulative_elements.back());
        } else {
            return sanisizer::create<std::vector<int> >(cumulative_elements.back());
        }
    }();

    for (I<decltype(num_barcodes)> i = 0; i < num_barcodes; ++i) {
        if (barcodes[i] > 0 && features[i] > 0 && features[i] <= nrow) {
            auto& oo = offsets[barcodes[i] - 1];
            const auto feat_idx = features[i] - 1;
            if constexpr(use_counts) {
                re_features[oo].first = feat_idx;
                re_features[oo].second = counts[i];
            } else {
                re_features[oo] = feat_idx;
            }
            ++oo;
        }
    }

    // Now counting the number of unique non-zero elements.
    auto nnz = sanisizer::create<std::vector<RLen> >(ncol);
    subpar::parallelize_range(num_threads, ncol, [&](int, int start, int length) {
        for (int c = start, end = start + length; c < end; ++c) {
            const auto first = cumulative_elements[c], last = cumulative_elements[c + 1];
            if (first == last) {
                continue;
            }

            std::sort(re_features.begin() + first, re_features.begin() + last);

            RLen count = 1;
            for (RLen i = first + 1; i < last; ++i) {
                const auto& current = re_features[i];
                const auto& previous = re_features[i - 1];
                if constexpr(use_counts) {
                    count += (current.first != previous.first);
                } else {
                    count += (current != previous);
                }
            }
            nnz[c] = count;
        }
    });

    auto populate_columns = [&](const auto& iptrs, const auto& vptrs) -> void {
        subpar::parallelize_range(num_threads, ncol, [&](int, int start, int length) {
            for (int c = start, end = start + length; c < end; ++c) {
                const auto first = cumulative_elements[c], last = cumulative_elements[c + 1];
                if (first == last) {
                    continue;
                }

                const auto iptr = iptrs[c];
                const auto vptr = vptrs[c];

                const auto& first_feature = re_features[first];
                if constexpr(use_counts) {
                    iptr[0] = first_feature.first;
                    vptr[0] = first_feature.second;
                } else {
                    sanisizer::cast<int>(nnz[c]);
                    iptr[0] = first_feature;
                    vptr[0] = 1;
                }

                RLen pos = 0;
                for (RLen i = first + 1; i < last; ++i) {
                    const auto& current = re_features[i];
                    const auto& previous = re_features[i - 1];
                    if constexpr(use_counts) {
                        if (current.first != previous.first) {
                            ++pos;
                            iptr[pos] = current.first;
                            vptr[pos] = current.second;
                        } else if constexpr(std::is_integral<I<decltype(vptr[pos])> >::value) {
                            vptr[pos] = sanisizer::sum<int>(vptr[pos], current.second);
                        } else {
                            vptr[pos] += current.second;
                        }
                    } else {
                        if (current != previous) {
                            ++pos;
                            iptr[pos] = current;
                            vptr[pos] = 1;
                        } else {
                            vptr[pos] += 1; // this is known to be safe, see cast check above.
                        }
                    }
                }
            }
        });
    };

    if (create_svt) {
        std::vector<Rcpp::IntegerVector> indices, values;
        indices.reserve(ncol);
        values.reserve(ncol);
        for (int c = 0; c < ncol; ++c) {
            indices.emplace_back(sanisizer::cast<I<decltype(std::declval<Rcpp::IntegerVector>().size())> >(nnz[c]));
            values.emplace_back(sanisizer::cast<I<decltype(std::declval<Rcpp::IntegerVector>().size())> >(nnz[c]));
        }

        std::vector<int*> iptrs, vptrs; 
        iptrs.reserve(ncol);
        vptrs.reserve(ncol);
        for (int c = 0; c < ncol; ++c) {
            iptrs.push_back(indices[c].begin());
            vptrs.push_back(values[c].begin());
        }
        populate_columns(iptrs, vptrs);

        auto output = sanisizer::create<Rcpp::List>(ncol);
        for (int c = 0; c < ncol; ++c) {
            output[c] = Rcpp::List::create(values[c], indices[c]);
        }
        return output;

    } else {
        std::vector<RLen> cumulative_elements2(sanisizer::sum<typename std::vector<RLen>::size_type>(ncol, 1));
        for (int c = 1; c <= ncol; ++c) {
            cumulative_elements2[c] = cumulative_elements2[c - 1] + nnz[c - 1];
        }

        auto indices = sanisizer::create<Rcpp::IntegerVector>(cumulative_elements2.back());
        auto values = sanisizer::create<Rcpp::NumericVector>(cumulative_elements2.back());

        std::vector<int*> iptrs;
        iptrs.reserve(ncol);
        std::vector<double*> vptrs; 
        vptrs.reserve(ncol);
        for (int c = 0; c < ncol; ++c) {
            iptrs.push_back(indices.begin() + cumulative_elements2[c]);
            vptrs.push_back(values.begin() + cumulative_elements2[c]);
        }
        populate_columns(iptrs, vptrs);

        auto ptrs = sanisizer::create<Rcpp::IntegerVector>(cumulative_elements2.size());
        for (int c = 1; c <= ncol; ++c) {
            ptrs[c] = sanisizer::cast<int>(cumulative_elements2[c]);
        }

        return Rcpp::List::create(std::move(indices), std::move(values), std::move(ptrs));
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::RObject make_count_matrix(
    int nrow,
    int ncol,
    Rcpp::IntegerVector features,
    Rcpp::IntegerVector barcodes,
    Rcpp::Nullable<Rcpp::IntegerVector> counts,
    bool create_svt,
    int num_threads
) {
    if (counts.isNull()) {
        return make_count_matrix_raw(nrow, ncol, std::move(features), std::move(barcodes), false, create_svt, num_threads);
    } else {
        return make_count_matrix_raw(nrow, ncol, std::move(features), std::move(barcodes), Rcpp::IntegerVector(counts), create_svt, num_threads);
    }
}
