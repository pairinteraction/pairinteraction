#include "pintr/diagonalizer/diagonalize.hpp"

#include "pintr/system/SystemAtom.hpp"

#include <complex>
#include <oneapi/tbb.h>

namespace pintr {
template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision) {
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range(systems.begin(), systems.end()),
                              [&](const auto &range) {
                                  for (auto &system : range) {
                                      system.get().diagonalize(diagonalizer, precision);
                                  }
                              });
}

template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 typename Derived::real_t min_eigenvalue, typename Derived::real_t max_eigenvalue,
                 int precision) {
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range(systems.begin(), systems.end()), [&](const auto &range) {
            for (auto &system : range) {
                system.get().diagonalize(diagonalizer, min_eigenvalue, max_eigenvalue, precision);
            }
        });
}

template <typename Derived>
void diagonalize(std::vector<Derived> &systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision) {
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range(systems.begin(), systems.end()),
                              [&](const auto &range) {
                                  for (auto &system : range) {
                                      system.diagonalize(diagonalizer, precision);
                                  }
                              });
}

template <typename Derived>
void diagonalize(std::vector<Derived> &systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 typename Derived::real_t min_eigenvalue, typename Derived::real_t max_eigenvalue,
                 int precision) {
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range(systems.begin(), systems.end()), [&](const auto &range) {
            for (auto &system : range) {
                system.diagonalize(diagonalizer, min_eigenvalue, max_eigenvalue, precision);
            }
        });
}

// Explicit instantiation
template void diagonalize(std::initializer_list<std::reference_wrapper<SystemAtom<float>>> systems,
                          const DiagonalizerInterface<SystemAtom<float>::scalar_t> &diagonalizer,
                          int precision);
template void diagonalize(std::initializer_list<std::reference_wrapper<SystemAtom<double>>> systems,
                          const DiagonalizerInterface<SystemAtom<double>::scalar_t> &diagonalizer,
                          int precision);
template void
diagonalize(std::initializer_list<std::reference_wrapper<SystemAtom<std::complex<float>>>> systems,
            const DiagonalizerInterface<SystemAtom<std::complex<float>>::scalar_t> &diagonalizer,
            int precision);
template void
diagonalize(std::initializer_list<std::reference_wrapper<SystemAtom<std::complex<double>>>> systems,
            const DiagonalizerInterface<SystemAtom<std::complex<double>>::scalar_t> &diagonalizer,
            int precision);

template void diagonalize(std::initializer_list<std::reference_wrapper<SystemAtom<float>>> systems,
                          const DiagonalizerInterface<SystemAtom<float>::scalar_t> &diagonalizer,
                          SystemAtom<float>::real_t min_eigenvalue,
                          SystemAtom<float>::real_t max_eigenvalue, int precision);
template void diagonalize(std::initializer_list<std::reference_wrapper<SystemAtom<double>>> systems,
                          const DiagonalizerInterface<SystemAtom<double>::scalar_t> &diagonalizer,
                          SystemAtom<double>::real_t min_eigenvalue,
                          SystemAtom<double>::real_t max_eigenvalue, int precision);
template void
diagonalize(std::initializer_list<std::reference_wrapper<SystemAtom<std::complex<float>>>> systems,
            const DiagonalizerInterface<SystemAtom<std::complex<float>>::scalar_t> &diagonalizer,
            SystemAtom<std::complex<float>>::real_t min_eigenvalue,
            SystemAtom<std::complex<float>>::real_t max_eigenvalue, int precision);
template void
diagonalize(std::initializer_list<std::reference_wrapper<SystemAtom<std::complex<double>>>> systems,
            const DiagonalizerInterface<SystemAtom<std::complex<double>>::scalar_t> &diagonalizer,
            SystemAtom<std::complex<double>>::real_t min_eigenvalue,
            SystemAtom<std::complex<double>>::real_t max_eigenvalue, int precision);

template void diagonalize(std::vector<SystemAtom<float>> &systems,
                          const DiagonalizerInterface<SystemAtom<float>::scalar_t> &diagonalizer,
                          int precision);
template void diagonalize(std::vector<SystemAtom<double>> &systems,
                          const DiagonalizerInterface<SystemAtom<double>::scalar_t> &diagonalizer,
                          int precision);
template void
diagonalize(std::vector<SystemAtom<std::complex<float>>> &systems,
            const DiagonalizerInterface<SystemAtom<std::complex<float>>::scalar_t> &diagonalizer,
            int precision);
template void
diagonalize(std::vector<SystemAtom<std::complex<double>>> &systems,
            const DiagonalizerInterface<SystemAtom<std::complex<double>>::scalar_t> &diagonalizer,
            int precision);

template void diagonalize(std::vector<SystemAtom<float>> &systems,
                          const DiagonalizerInterface<SystemAtom<float>::scalar_t> &diagonalizer,
                          SystemAtom<float>::real_t min_eigenvalue,
                          SystemAtom<float>::real_t max_eigenvalue, int precision);
template void diagonalize(std::vector<SystemAtom<double>> &systems,
                          const DiagonalizerInterface<SystemAtom<double>::scalar_t> &diagonalizer,
                          SystemAtom<double>::real_t min_eigenvalue,
                          SystemAtom<double>::real_t max_eigenvalue, int precision);
template void
diagonalize(std::vector<SystemAtom<std::complex<float>>> &systems,
            const DiagonalizerInterface<SystemAtom<std::complex<float>>::scalar_t> &diagonalizer,
            SystemAtom<std::complex<float>>::real_t min_eigenvalue,
            SystemAtom<std::complex<float>>::real_t max_eigenvalue, int precision);
template void
diagonalize(std::vector<SystemAtom<std::complex<double>>> &systems,
            const DiagonalizerInterface<SystemAtom<std::complex<double>>::scalar_t> &diagonalizer,
            SystemAtom<std::complex<double>>::real_t min_eigenvalue,
            SystemAtom<std::complex<double>>::real_t max_eigenvalue, int precision);
} // namespace pintr
