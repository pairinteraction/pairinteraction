// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/diagonalize/diagonalize.hpp"

#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/system/SystemPair.hpp"
#include "pairinteraction/utils/Range.hpp"

#include <complex>
#include <oneapi/tbb.h>
#include <optional>

namespace pairinteraction {

template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 std::optional<typename Derived::real_t> min_eigenenergy,
                 std::optional<typename Derived::real_t> max_eigenenergy, double rtol) {
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range(systems.begin(), systems.end()), [&](const auto &range) {
            for (auto &system : range) {
                system.get().diagonalize(diagonalizer, min_eigenenergy, max_eigenenergy, rtol);
            }
        });
}

template <typename Derived>
void diagonalize(std::vector<Derived> &systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 std::optional<typename Derived::real_t> min_eigenenergy,
                 std::optional<typename Derived::real_t> max_eigenenergy, double rtol) {
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range(systems.begin(), systems.end()), [&](const auto &range) {
            for (auto &system : range) {
                system.diagonalize(diagonalizer, min_eigenenergy, max_eigenenergy, rtol);
            }
        });
}

template <typename Derived>
void diagonalize(std::vector<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 std::optional<typename Derived::real_t> min_eigenenergy,
                 std::optional<typename Derived::real_t> max_eigenenergy, double rtol) {
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range(systems.begin(), systems.end()), [&](const auto &range) {
            for (auto &system : range) {
                system.get().diagonalize(diagonalizer, min_eigenenergy, max_eigenenergy, rtol);
            }
        });
}

// Explicit instantiations
// NOLINTBEGIN(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)
#define INSTANTIATE_DIAGONALIZE_HELPER(SCALAR, TYPE)                                               \
    template void diagonalize(std::initializer_list<std::reference_wrapper<TYPE<SCALAR>>> systems, \
                              const DiagonalizerInterface<TYPE<SCALAR>::scalar_t> &diagonalizer,   \
                              std::optional<TYPE<SCALAR>::real_t> min_eigenenergy,                 \
                              std::optional<TYPE<SCALAR>::real_t> max_eigenenergy, double rtol);   \
    template void diagonalize(std::vector<TYPE<SCALAR>> &systems,                                  \
                              const DiagonalizerInterface<TYPE<SCALAR>::scalar_t> &diagonalizer,   \
                              std::optional<TYPE<SCALAR>::real_t> min_eigenenergy,                 \
                              std::optional<TYPE<SCALAR>::real_t> max_eigenenergy, double rtol);   \
    template void diagonalize(std::vector<std::reference_wrapper<TYPE<SCALAR>>> systems,           \
                              const DiagonalizerInterface<TYPE<SCALAR>::scalar_t> &diagonalizer,   \
                              std::optional<TYPE<SCALAR>::real_t> min_eigenenergy,                 \
                              std::optional<TYPE<SCALAR>::real_t> max_eigenenergy, double rtol);
#define INSTANTIATE_DIAGONALIZE(SCALAR)                                                            \
    INSTANTIATE_DIAGONALIZE_HELPER(SCALAR, SystemAtom)                                             \
    INSTANTIATE_DIAGONALIZE_HELPER(SCALAR, SystemPair)
// NOLINTEND(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)

INSTANTIATE_DIAGONALIZE(double)
INSTANTIATE_DIAGONALIZE(std::complex<double>)

#undef INSTANTIATE_DIAGONALIZE_HELPER
#undef INSTANTIATE_DIAGONALIZE

} // namespace pairinteraction
