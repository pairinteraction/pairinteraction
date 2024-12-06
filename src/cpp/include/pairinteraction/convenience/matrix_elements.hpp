#pragma once

#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <memory>

namespace pairinteraction {

template <typename Real>
class KetAtom;

template <typename Scalar>
class SystemAtom;

enum class OperatorType;

template <typename Scalar>
typename traits::NumTraits<Scalar>::real_t
calculate_energy(std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> ket,
                 const SystemAtom<Scalar> &system);

template <typename Real>
Real calculate_energy(std::shared_ptr<const KetAtom<Real>> ket);

template <typename Scalar>
Scalar calculate_matrix_element(
    std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> initial_ket,
    std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> final_ket,
    const SystemAtom<Scalar> &system, OperatorType type, int q);

template <typename Real>
Real calculate_matrix_element(std::shared_ptr<const KetAtom<Real>> initial_ket,
                              std::shared_ptr<const KetAtom<Real>> final_ket, OperatorType type,
                              int q);

template <typename Scalar>
Scalar calculate_electric_dipole_matrix_element(
    std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> initial_ket,
    std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> final_ket,
    const SystemAtom<Scalar> &system, int q);

template <typename Real>
Real calculate_electric_dipole_matrix_element(std::shared_ptr<const KetAtom<Real>> initial_ket,
                                              std::shared_ptr<const KetAtom<Real>> final_ket,
                                              int q);

// Extern template declarations
// NOLINTBEGIN(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)
#define EXTERN_CALCULATE_FUNCTIONS_REAL(REAL)                                                      \
    extern template REAL calculate_energy<REAL>(std::shared_ptr<const KetAtom<REAL>>);             \
    extern template REAL calculate_matrix_element<REAL>(std::shared_ptr<const KetAtom<REAL>>,      \
                                                        std::shared_ptr<const KetAtom<REAL>>,      \
                                                        OperatorType, int);                        \
    extern template REAL calculate_electric_dipole_matrix_element<REAL>(                           \
        std::shared_ptr<const KetAtom<REAL>>, std::shared_ptr<const KetAtom<REAL>>, int);

#define EXTERN_CALCULATE_FUNCTIONS(SCALAR, REAL)                                                   \
    extern template REAL calculate_energy<SCALAR>(std::shared_ptr<const KetAtom<REAL>>,            \
                                                  const SystemAtom<SCALAR> &);                     \
    extern template SCALAR calculate_matrix_element<SCALAR>(                                       \
        std::shared_ptr<const KetAtom<REAL>>, std::shared_ptr<const KetAtom<REAL>>,                \
        const SystemAtom<SCALAR> &, OperatorType, int);                                            \
    extern template SCALAR calculate_electric_dipole_matrix_element<SCALAR>(                       \
        std::shared_ptr<const KetAtom<REAL>>, std::shared_ptr<const KetAtom<REAL>>,                \
        const SystemAtom<SCALAR> &, int);
// NOLINTEND(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)

EXTERN_CALCULATE_FUNCTIONS_REAL(float)
EXTERN_CALCULATE_FUNCTIONS_REAL(double)

EXTERN_CALCULATE_FUNCTIONS(float, float)
EXTERN_CALCULATE_FUNCTIONS(double, double)
EXTERN_CALCULATE_FUNCTIONS(std::complex<float>, float)
EXTERN_CALCULATE_FUNCTIONS(std::complex<double>, double)

#undef EXTERN_CALCULATE_FUNCTIONS_REAL
#undef EXTERN_CALCULATE_FUNCTIONS
} // namespace pairinteraction
