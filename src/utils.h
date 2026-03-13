#ifndef UTILS_H
#define UTILS_H

#include <type_traits>

template <typename Type_>
using I = std::remove_cv_t<std::remove_reference_t<Type_> >;

#endif
