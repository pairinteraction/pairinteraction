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
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/SparseCore>
#include <memory>

namespace pairinteraction {
template <typename Scalar>
SystemCombined<Scalar>::SystemCombined(const SystemAtom<Scalar> &system1,
                                       const SystemAtom<Scalar> &system2, real_t min_energy,
                                       real_t max_energy)
    : System<SystemCombined<Scalar>>(
          create_combined_basis(system1, system2, min_energy, max_energy)),
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
        auto d1_plus = OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, 1);
        auto d1_minus = OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, -1);
        auto d1_zero = OperatorAtom<Scalar>(basis1, OperatorType::ELECTRIC_DIPOLE, 0);
        auto d2_plus = OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, 1);
        auto d2_minus = OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, -1);
        auto d2_zero = OperatorAtom<Scalar>(basis2, OperatorType::ELECTRIC_DIPOLE, 0);

        // Test: combine two operator matrices
        real_t precision = 10 * std::numeric_limits<real_t>::epsilon();

        const auto &matrix1 = d1_plus.get_matrix();
        const auto &matrix2 = d2_minus.get_matrix();

        std::vector<Eigen::Triplet<Scalar>> triplets;

        // Loop over the rows of the first matrix
        for (Eigen::Index outer_idx1 = 0; outer_idx1 < matrix1.outerSize(); ++outer_idx1) {

            // Loop over the non-zero column elements of the first matrix
            for (typename Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator it1(
                     matrix1, outer_idx1);
                 it1; ++it1) {

                Eigen::Index row1 = it1.row();
                Eigen::Index col1 = it1.col();
                Scalar value1 = it1.value();

                const auto &range_outer_idx2 =
                    this->hamiltonian->get_basis()->get_index_range(outer_idx1);
                const auto &range_inner_idx2 =
                    this->hamiltonian->get_basis()->get_index_range(it1.index());

                // Loop over the rows of the second matrix that are energetically allowed
                for (Eigen::Index outer_idx2 = static_cast<Eigen::Index>(range_outer_idx2.min());
                     outer_idx2 < static_cast<Eigen::Index>(range_outer_idx2.max()); ++outer_idx2) {

                    // Calculate the minimum and maximum values of the index pointer of the second
                    // matrix
                    Eigen::Index begin_idxptr2 = matrix2.outerIndexPtr()[outer_idx2];
                    Eigen::Index end_idxptr2 = matrix2.outerIndexPtr()[outer_idx2 + 1];

                    // The minimum value is chosen such that we start with an energetically allowed
                    // column
                    begin_idxptr2 +=
                        std::distance(matrix2.innerIndexPtr() + begin_idxptr2,
                                      std::lower_bound(matrix2.innerIndexPtr() + begin_idxptr2,
                                                       matrix2.innerIndexPtr() + end_idxptr2,
                                                       range_inner_idx2.min()));

                    // Loop over the non-zero column elements of the second matrix that are
                    // energetically allowed (we break the loop if the index pointer corresponds to
                    // a column that is not energetically allowed)
                    for (Eigen::Index idxptr2 = begin_idxptr2; idxptr2 < end_idxptr2; ++idxptr2) {

                        Eigen::Index col2 = matrix2.innerIndexPtr()[idxptr2];
                        if (col2 >= static_cast<Eigen::Index>(range_inner_idx2.max())) {
                            break;
                        }
                        Eigen::Index row2 = outer_idx2;
                        Scalar value2 = matrix2.valuePtr()[idxptr2];

                        // Calculate the row and column index of the entry in the combined matrix
                        Eigen::Index row =
                            this->hamiltonian->get_basis()->get_combined_index(row1, row2);
                        Eigen::Index col =
                            this->hamiltonian->get_basis()->get_combined_index(col1, col2);

                        // Calculate the value of the entry in the combined matrix
                        Scalar value = value1 * value2;

                        // Store the entry
                        if (std::abs(value) > precision) {
                            triplets.emplace_back(row, col, value);
                        }
                    }
                }
            }
        }

        // Construct the combined matrix from the triplets
        Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix(basis->get_number_of_states(),
                                                            basis->get_number_of_states());
        matrix.setFromTriplets(triplets.begin(), triplets.end());
        matrix.makeCompressed();

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
std::shared_ptr<const typename SystemCombined<Scalar>::basis_t>
SystemCombined<Scalar>::create_combined_basis(const SystemAtom<Scalar> &system1,
                                              const SystemAtom<Scalar> &system2, real_t min_energy,
                                              real_t max_energy) const {

    // Get the eigenvalues of the constituent systems
    auto eigenvalues1 = system1.get_eigenvalues();
    auto eigenvalues2 = system2.get_eigenvalues();

    // Construct the canonical basis that contains all all combined states with energies between
    // min_energy and max_energy
    ketvec_t kets;
    kets.reserve(eigenvalues1.size() * eigenvalues2.size());

    typename basis_t::map_size_t map_index_combined_state;
    map_index_combined_state.reserve(eigenvalues1.size() * eigenvalues2.size());

    typename basis_t::map_range_t map_range_of_index_state2;
    map_range_of_index_state2.reserve(eigenvalues1.size());

    // Loop only over states with an allowed energy // TODO make this more efficient
    for (size_t idx1 = 0; idx1 < static_cast<size_t>(eigenvalues1.size()); ++idx1) {

        // TODO set the following properly (requires the states to be sorted by energy)
        map_range_of_index_state2.emplace(idx1, typename basis_t::range_t(0, eigenvalues2.size()));

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

            // Store the combined index
            map_index_combined_state[idx1 * eigenvalues2.size() + idx2] = kets.size();

            // Store the combined state as a ket
            kets.emplace_back(std::make_shared<ket_t>(
                kets.size(), energy, quantum_number_f, quantum_number_m, parity,
                std::vector<std::shared_ptr<const KetAtom<real_t>>>{ket1, ket2}));
        }
    }

    kets.shrink_to_fit();

    return std::make_shared<BasisCombined<Scalar>>(
        std::move(kets), std::move(map_index_combined_state), std::move(map_range_of_index_state2),
        eigenvalues2.size());
}

// Explicit instantiations
template class SystemCombined<float>;
template class SystemCombined<double>;
template class SystemCombined<std::complex<float>>;
template class SystemCombined<std::complex<double>>;
} // namespace pairinteraction
