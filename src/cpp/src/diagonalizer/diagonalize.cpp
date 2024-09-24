#include "pairinteraction/diagonalizer/diagonalize.hpp"

#include "pairinteraction/system/SystemAtom.hpp"
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
// NOLINTNEXTLINE(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)
#define INSTANTIATE_DIAGONALIZE(TYPE)                                                              \
    template void diagonalize(                                                                     \
        std::initializer_list<std::reference_wrapper<SystemAtom<TYPE>>> systems,                   \
        const DiagonalizerInterface<SystemAtom<TYPE>::scalar_t> &diagonalizer, int precision,      \
        const Range<SystemAtom<TYPE>::real_t> &eigenvalue_range);                                  \
    template void diagonalize(                                                                     \
        std::vector<SystemAtom<TYPE>> &systems,                                                    \
        const DiagonalizerInterface<SystemAtom<TYPE>::scalar_t> &diagonalizer, int precision,      \
        const Range<SystemAtom<TYPE>::real_t> &eigenvalue_range);                                  \
    template void diagonalize(                                                                     \
        std::vector<std::reference_wrapper<SystemAtom<TYPE>>> systems,                             \
        const DiagonalizerInterface<SystemAtom<TYPE>::scalar_t> &diagonalizer, int precision,      \
        const Range<SystemAtom<TYPE>::real_t> &eigenvalue_range);

INSTANTIATE_DIAGONALIZE(float)
INSTANTIATE_DIAGONALIZE(double)
INSTANTIATE_DIAGONALIZE(std::complex<float>)
INSTANTIATE_DIAGONALIZE(std::complex<double>)

#undef INSTANTIATE_DIAGONALIZE

} // namespace pairinteraction
