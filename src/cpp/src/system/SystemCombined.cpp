#include "pairinteraction/system/SystemCombined.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisCombined.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetCombined.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/operator/OperatorCombined.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/SparseCore>
#include <memory>

namespace pairinteraction {
template <typename Scalar>
SystemCombined<Scalar>::SystemCombined(const SystemAtom<Scalar> &system1,
                                       const SystemAtom<Scalar> &system2, real_t min_energy,
                                       real_t max_energy)
    : System<SystemCombined<Scalar>>(std::make_shared<BasisCombined<Scalar>>(
          create_combined_kets(system1, system2, min_energy, max_energy))),
      basis1(system1.get_basis()), basis2(system2.get_basis()) {}

template <typename Scalar>
SystemCombined<Scalar> &SystemCombined<Scalar>::set_interatomic_distance(real_t distance) {
    this->interatomic_distance = distance;
    this->hamiltonian_requires_construction = true;
    return *this;
}

template <typename Scalar>
void SystemCombined<Scalar>::construct_hamiltonian() const {
    auto basis = this->hamiltonian->get_basis();

    // Construct the unperturbed Hamiltonian
    this->hamiltonian = std::make_unique<OperatorCombined<Scalar>>(basis, OperatorType::ENERGY);
    this->hamiltonian_is_diagonal = true;

    // TODO if single-atom stats can be labeled by the quantum number m, the following line
    // should set the variable to true
    bool sort_by_quantum_number_m = false;

    // Dipole-dipole interaction along the z-axis
    if (this->interatomic_distance != 0) {

        // Get the constituting single-atom operators
        auto d1_plus = OperatorAtom<Scalar>(basis1, OperatorType::MAGNETIC_DIPOLE, 1);
        auto d1_minus = OperatorAtom<Scalar>(basis1, OperatorType::MAGNETIC_DIPOLE, -1);
        auto d1_zero = OperatorAtom<Scalar>(basis1, OperatorType::MAGNETIC_DIPOLE, 0);
        auto d2_plus = OperatorAtom<Scalar>(basis2, OperatorType::MAGNETIC_DIPOLE, 1);
        auto d2_minus = OperatorAtom<Scalar>(basis2, OperatorType::MAGNETIC_DIPOLE, -1);
        auto d2_zero = OperatorAtom<Scalar>(basis2, OperatorType::MAGNETIC_DIPOLE, 0);

        // TODO build the dipole-dipole interaction

        this->hamiltonian_is_diagonal = false;
    }

    // Store which labels can be used to block-diagonalize the Hamiltonian
    this->blockdiagonalizing_labels.clear();
    if (sort_by_quantum_number_m) {
        this->blockdiagonalizing_labels.push_back(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    }
}

template <typename Scalar>
typename SystemCombined<Scalar>::ketvec_t
SystemCombined<Scalar>::create_combined_kets(const SystemAtom<Scalar> &system1,
                                             const SystemAtom<Scalar> &system2, real_t min_energy,
                                             real_t max_energy) const {
    // Get the eigenvalues of the constituent systems
    auto eigenvalues1 = system1.get_eigenvalues();
    auto eigenvalues2 = system2.get_eigenvalues();

    // Construct the canonical basis that contains all all combined states with energies between
    // min_energy and max_energy
    ketvec_t kets;
    kets.reserve(eigenvalues1.size() * eigenvalues2.size());

    // Loop only over states with an allowed energy // TODO make this more efficient
    size_t id = 0;
    for (size_t idx1 = 0; idx1 < static_cast<size_t>(eigenvalues1.size()); ++idx1) {
        for (size_t idx2 = 0; idx2 < static_cast<size_t>(eigenvalues2.size()); ++idx2) {
            const real_t energy = eigenvalues1[idx1] + eigenvalues2[idx2];
            if (energy < min_energy || energy > max_energy) {
                continue;
            }

            // Get quantum numbers
            // TODO symmetrize the states so that parity and quantum_number_f are also well-defined
            Parity parity = Parity::UNKNOWN;
            real_t quantum_number_f = std::numeric_limits<real_t>::max();
            real_t quantum_number_m = std::numeric_limits<real_t>::max();
            try {
                // TODO make it work without a try-catch block even if the quantum numbers are not
                // defined
                quantum_number_m = system1.get_basis()->get_quantum_number_m(idx1) +
                    system2.get_basis()->get_quantum_number_m(idx2);
            } catch (const std::invalid_argument &e) {
            }

            // Get ket with largest overlap
            std::shared_ptr<const KetAtom<real_t>> ket1 =
                system1.get_basis()->get_ket_with_largest_overlap(idx1);
            std::shared_ptr<const KetAtom<real_t>> ket2 =
                system2.get_basis()->get_ket_with_largest_overlap(idx2);

            // Store the combined state as a ket
            kets.emplace_back(std::make_shared<ket_t>(
                id, energy, quantum_number_f, quantum_number_m, parity,
                std::vector<std::shared_ptr<const KetAtom<real_t>>>{ket1, ket2}));
            ++id;
        }
    }

    kets.shrink_to_fit();
    return kets;
}

// Explicit instantiations
template class SystemCombined<float>;
template class SystemCombined<double>;
template class SystemCombined<std::complex<float>>;
template class SystemCombined<std::complex<double>>;
} // namespace pairinteraction
