#include "pintr/interfaces/DiagonalizerInterface.hpp"

#include <thread>

namespace pintr {
template <typename Scalar>
DiagonalizerInterface<Scalar>::DiagonalizerInterface(int num_cpu_cores)
    : smph_cpu_cores(num_cpu_cores > 0 ? num_cpu_cores
                                       : std::max(std::thread::hardware_concurrency(), 1U)) {}

template class DiagonalizerInterface<float>;
template class DiagonalizerInterface<double>;
template class DiagonalizerInterface<std::complex<float>>;
template class DiagonalizerInterface<std::complex<double>>;
} // namespace pintr
