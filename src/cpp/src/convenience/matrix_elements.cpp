#include "pairinteraction/convenience/matrix_elements.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <memory>

namespace pairinteraction {

template <typename Scalar>
typename traits::NumTraits<Scalar>::real_t
calculate_energy(std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> ket,
                 const SystemAtom<Scalar> &system) {

    auto state = system.get_basis()->get_state_with_largest_overlap(ket);

    // TODO throw an exception if ket has a small overlap with the state

    auto matrix_element =
        ket->get_database().get_matrix_elements(state, state, OperatorType::ENERGY, 0);

    assert(matrix_element.rows() == 1);
    assert(matrix_element.cols() == 1);

    return std::real(matrix_element.coeff(0, 0));
}

template <typename Real>
Real calculate_energy(std::shared_ptr<const KetAtom<Real>> ket) {
    return ket->get_energy();
}

template <typename Scalar>
Scalar calculate_electric_dipole_matrix_element(
    std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> initial_ket,
    std::shared_ptr<const KetAtom<typename traits::NumTraits<Scalar>::real_t>> final_ket,
    const SystemAtom<Scalar> &system, int q) {
    auto initial_state = system.get_basis()->get_state_with_largest_overlap(initial_ket);
    auto final_state = system.get_basis()->get_state_with_largest_overlap(final_ket);

    // TODO throw an exception if a ket has a small overlap with a state

    auto matrix_element = initial_ket->get_database().get_matrix_elements(
        initial_state, final_state, OperatorType::ELECTRIC_DIPOLE, q);

    assert(matrix_element.rows() == 1);
    assert(matrix_element.cols() == 1);

    return matrix_element.coeff(0, 0);
}

template <typename Real>
std::complex<Real>
calculate_electric_dipole_matrix_element(std::shared_ptr<const KetAtom<Real>> initial_ket,
                                         std::shared_ptr<const KetAtom<Real>> final_ket, int q) {
    auto basis = BasisAtomCreator<Real>()
                     .append_ket(initial_ket)
                     .append_ket(final_ket)
                     .create(initial_ket->get_database());
    SystemAtom<Real> system(basis);

    return calculate_electric_dipole_matrix_element(initial_ket, final_ket, system, q);
}

// TODO add test cases in cpp and python for this calculate functions

// Explicit instantiation
// TODO simplify the following code using macros
template float calculate_energy<float>(std::shared_ptr<const KetAtom<float>>,
                                       const SystemAtom<float> &);
template double calculate_energy<double>(std::shared_ptr<const KetAtom<double>>,
                                         const SystemAtom<double> &);
template float calculate_energy<std::complex<float>>(std::shared_ptr<const KetAtom<float>>,
                                                     const SystemAtom<std::complex<float>> &);
template double calculate_energy<std::complex<double>>(std::shared_ptr<const KetAtom<double>>,
                                                       const SystemAtom<std::complex<double>> &);

template float
calculate_electric_dipole_matrix_element<float>(std::shared_ptr<const KetAtom<float>>,
                                                std::shared_ptr<const KetAtom<float>>,
                                                const SystemAtom<float> &, int);
template double
calculate_electric_dipole_matrix_element<double>(std::shared_ptr<const KetAtom<double>>,
                                                 std::shared_ptr<const KetAtom<double>>,
                                                 const SystemAtom<double> &, int);
template std::complex<float> calculate_electric_dipole_matrix_element<std::complex<float>>(
    std::shared_ptr<const KetAtom<float>>, std::shared_ptr<const KetAtom<float>>,
    const SystemAtom<std::complex<float>> &, int);
template std::complex<double> calculate_electric_dipole_matrix_element<std::complex<double>>(
    std::shared_ptr<const KetAtom<double>>, std::shared_ptr<const KetAtom<double>>,
    const SystemAtom<std::complex<double>> &, int);

template float calculate_energy<float>(std::shared_ptr<const KetAtom<float>>);
template double calculate_energy<double>(std::shared_ptr<const KetAtom<double>>);

template std::complex<float>
calculate_electric_dipole_matrix_element<float>(std::shared_ptr<const KetAtom<float>>,
                                                std::shared_ptr<const KetAtom<float>>, int);
template std::complex<double>
calculate_electric_dipole_matrix_element<double>(std::shared_ptr<const KetAtom<double>>,
                                                 std::shared_ptr<const KetAtom<double>>, int);

} // namespace pairinteraction
