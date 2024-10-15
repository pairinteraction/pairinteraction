#include "pairinteraction/diagonalizer/diagonalize.hpp"

#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/system/SystemCombined.hpp"
#include "pairinteraction/utils/Range.hpp"

#include <complex>
#include <oneapi/tbb.h>

namespace pairinteraction {

template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision, const Range<typename Derived::real_t> &eigenvalue_range) {
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range(systems.begin(), systems.end()), [&](const auto &range) {
            for (auto &system : range) {
                system.get().diagonalize(diagonalizer, precision, eigenvalue_range);
            }
        });
}

template <typename Derived>
void diagonalize(std::vector<Derived> &systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision, const Range<typename Derived::real_t> &eigenvalue_range) {
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range(systems.begin(), systems.end()),
                              [&](const auto &range) {
                                  for (auto &system : range) {
                                      system.diagonalize(diagonalizer, precision, eigenvalue_range);
                                  }
                              });
}

template <typename Derived>
void diagonalize(std::vector<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision, const Range<typename Derived::real_t> &eigenvalue_range) {
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range(systems.begin(), systems.end()), [&](const auto &range) {
            for (auto &system : range) {
                system.get().diagonalize(diagonalizer, precision, eigenvalue_range);
            }
        });
}

// Explicit instantiations
// NOLINTBEGIN(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)
#define INSTANTIATE_DIAGONALIZE_HELPER(SCALAR, TYPE)                                               \
    template void diagonalize(std::initializer_list<std::reference_wrapper<TYPE<SCALAR>>> systems, \
                              const DiagonalizerInterface<TYPE<SCALAR>::scalar_t> &diagonalizer,   \
                              int precision, const Range<TYPE<SCALAR>::real_t> &eigenvalue_range); \
    template void diagonalize(std::vector<TYPE<SCALAR>> &systems,                                  \
                              const DiagonalizerInterface<TYPE<SCALAR>::scalar_t> &diagonalizer,   \
                              int precision, const Range<TYPE<SCALAR>::real_t> &eigenvalue_range); \
    template void diagonalize(std::vector<std::reference_wrapper<TYPE<SCALAR>>> systems,           \
                              const DiagonalizerInterface<TYPE<SCALAR>::scalar_t> &diagonalizer,   \
                              int precision, const Range<TYPE<SCALAR>::real_t> &eigenvalue_range);
#define INSTANTIATE_DIAGONALIZE(SCALAR)                                                            \
    INSTANTIATE_DIAGONALIZE_HELPER(SCALAR, SystemAtom)                                             \
    INSTANTIATE_DIAGONALIZE_HELPER(SCALAR, SystemCombined)
// NOLINTEND(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)

INSTANTIATE_DIAGONALIZE(float)
INSTANTIATE_DIAGONALIZE(double)
INSTANTIATE_DIAGONALIZE(std::complex<float>)
INSTANTIATE_DIAGONALIZE(std::complex<double>)

#undef INSTANTIATE_DIAGONALIZE_HELPER
#undef INSTANTIATE_DIAGONALIZE

} // namespace pairinteraction
