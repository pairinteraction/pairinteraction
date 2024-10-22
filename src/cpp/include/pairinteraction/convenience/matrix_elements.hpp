#pragma once

#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <memory>

namespace pairinteraction {

template <typename Real>
class KetAtom;

template <typename Scalar>
class SystemAtom;

template <typename Scalar>
typename traits::NumTraits<Scalar>::real_t
calculate_energy(std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> ket,
                 const SystemAtom<Scalar> &system);

template <typename Real>
Real calculate_energy(std::shared_ptr<const KetAtom<Real>> ket);

template <typename Scalar>
Scalar calculate_electric_dipole_matrix_element(
    std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> initial_ket,
    std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> final_ket,
    const SystemAtom<Scalar> &system, int q);

template <typename Real>
Real calculate_electric_dipole_matrix_element(std::shared_ptr<const KetAtom<Real>> initial_ket,
                                              std::shared_ptr<const KetAtom<Real>> final_ket,
                                              int q);

// TODO simplify the following code using macros
extern template float calculate_energy<float>(std::shared_ptr<const KetAtom<float>>,
                                              const SystemAtom<float> &);
extern template double calculate_energy<double>(std::shared_ptr<const KetAtom<double>>,
                                                const SystemAtom<double> &);
extern template float
calculate_energy<std::complex<float>>(std::shared_ptr<const KetAtom<float>>,
                                      const SystemAtom<std::complex<float>> &);
extern template double
calculate_energy<std::complex<double>>(std::shared_ptr<const KetAtom<double>>,
                                       const SystemAtom<std::complex<double>> &);

extern template float
calculate_electric_dipole_matrix_element<float>(std::shared_ptr<const KetAtom<float>>,
                                                std::shared_ptr<const KetAtom<float>>,
                                                const SystemAtom<float> &, int);
extern template double
calculate_electric_dipole_matrix_element<double>(std::shared_ptr<const KetAtom<double>>,
                                                 std::shared_ptr<const KetAtom<double>>,
                                                 const SystemAtom<double> &, int);
extern template std::complex<float> calculate_electric_dipole_matrix_element<std::complex<float>>(
    std::shared_ptr<const KetAtom<float>>, std::shared_ptr<const KetAtom<float>>,
    const SystemAtom<std::complex<float>> &, int);
extern template std::complex<double> calculate_electric_dipole_matrix_element<std::complex<double>>(
    std::shared_ptr<const KetAtom<double>>, std::shared_ptr<const KetAtom<double>>,
    const SystemAtom<std::complex<double>> &, int);

extern template float calculate_energy<float>(std::shared_ptr<const KetAtom<float>>);
extern template double calculate_energy<double>(std::shared_ptr<const KetAtom<double>>);

extern template float
calculate_electric_dipole_matrix_element<float>(std::shared_ptr<const KetAtom<float>>,
                                                std::shared_ptr<const KetAtom<float>>, int);
extern template double
calculate_electric_dipole_matrix_element<double>(std::shared_ptr<const KetAtom<double>>,
                                                 std::shared_ptr<const KetAtom<double>>, int);

} // namespace pairinteraction
