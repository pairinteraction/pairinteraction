#include "pintr/diagonalizer/diagonalize.hpp"

#include "pintr/interfaces/TransformationBuilderInterface.hpp"
#include "pintr/system/System.hpp"
#include "pintr/system/SystemAtom.hpp"

#include <complex>
#include <oneapi/tbb.h>

namespace pintr {
template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<System<Derived>>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision) {
    using real_t = typename Derived::real_t;
    real_t min_eigenvalue = std::numeric_limits<real_t>::min();
    real_t max_eigenvalue = std::numeric_limits<real_t>::max();
    diagonalize(systems, diagonalizer, min_eigenvalue, max_eigenvalue, precision);
}

template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<System<Derived>>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 typename Derived::real_t min_eigenvalue, typename Derived::real_t max_eigenvalue,
                 int precision) {
    oneapi::tbb::parallel_for_each(systems.begin(), systems.end(), [&](System<Derived> &system) {
        system.diagonalize(diagonalizer, min_eigenvalue, max_eigenvalue, precision);
    });
}

template void
diagonalize(std::initializer_list<std::reference_wrapper<System<SystemAtom<float>>>> systems,
            const DiagonalizerInterface<SystemAtom<float>::scalar_t> &diagonalizer, int precision);
template void
diagonalize(std::initializer_list<std::reference_wrapper<System<SystemAtom<double>>>> systems,
            const DiagonalizerInterface<SystemAtom<double>::scalar_t> &diagonalizer, int precision);
template void diagonalize(
    std::initializer_list<std::reference_wrapper<System<SystemAtom<std::complex<float>>>>> systems,
    const DiagonalizerInterface<SystemAtom<std::complex<float>>::scalar_t> &diagonalizer,
    int precision);
template void diagonalize(
    std::initializer_list<std::reference_wrapper<System<SystemAtom<std::complex<double>>>>> systems,
    const DiagonalizerInterface<SystemAtom<std::complex<double>>::scalar_t> &diagonalizer,
    int precision);

template void
diagonalize(std::initializer_list<std::reference_wrapper<System<SystemAtom<float>>>> systems,
            const DiagonalizerInterface<SystemAtom<float>::scalar_t> &diagonalizer,
            SystemAtom<float>::real_t min_eigenvalue, SystemAtom<float>::real_t max_eigenvalue,
            int precision);
template void
diagonalize(std::initializer_list<std::reference_wrapper<System<SystemAtom<double>>>> systems,
            const DiagonalizerInterface<SystemAtom<double>::scalar_t> &diagonalizer,
            SystemAtom<double>::real_t min_eigenvalue, SystemAtom<double>::real_t max_eigenvalue,
            int precision);
template void diagonalize(
    std::initializer_list<std::reference_wrapper<System<SystemAtom<std::complex<float>>>>> systems,
    const DiagonalizerInterface<SystemAtom<std::complex<float>>::scalar_t> &diagonalizer,
    SystemAtom<std::complex<float>>::real_t min_eigenvalue,
    SystemAtom<std::complex<float>>::real_t max_eigenvalue, int precision);
template void diagonalize(
    std::initializer_list<std::reference_wrapper<System<SystemAtom<std::complex<double>>>>> systems,
    const DiagonalizerInterface<SystemAtom<std::complex<double>>::scalar_t> &diagonalizer,
    SystemAtom<std::complex<double>>::real_t min_eigenvalue,
    SystemAtom<std::complex<double>>::real_t max_eigenvalue, int precision);
} // namespace pintr
