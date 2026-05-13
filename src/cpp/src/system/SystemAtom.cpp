// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/SystemAtom.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/system/GreenTensorInterpolator.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/operator.hpp"
#include "pairinteraction/utils/spherical.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <memory>
#include <set>
#include <spdlog/spdlog.h>
#include <variant>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
SystemAtom<Scalar>::SystemAtom(std::shared_ptr<const basis_t> basis)
    : System<SystemAtom<Scalar>>(std::move(basis)) {}

template <typename Scalar>
SystemAtom<Scalar> &SystemAtom<Scalar>::set_electric_field(const std::array<real_t, 3> &field) {
    this->hamiltonian_requires_construction = true;

    if (!traits::NumTraits<Scalar>::is_complex_v && field[1] != 0) {
        throw std::invalid_argument(
            "The field must not have a y-component if the scalar type is real.");
    }

    electric_field = field;

    return *this;
}

template <typename Scalar>
SystemAtom<Scalar> &SystemAtom<Scalar>::set_magnetic_field(const std::array<real_t, 3> &field) {
    this->hamiltonian_requires_construction = true;

    if (!traits::NumTraits<Scalar>::is_complex_v && field[1] != 0) {
        throw std::invalid_argument(
            "The field must not have a y-component if the scalar type is real.");
    }

    magnetic_field = field;

    return *this;
}

template <typename Scalar>
SystemAtom<Scalar> &SystemAtom<Scalar>::set_diamagnetism_enabled(bool enable) {
    this->hamiltonian_requires_construction = true;
    diamagnetism_enabled = enable;
    return *this;
}

template <typename Scalar>
SystemAtom<Scalar> &
SystemAtom<Scalar>::set_ion_distance_vector(const std::array<real_t, 3> &vector) {
    this->hamiltonian_requires_construction = true;

    if (!traits::NumTraits<Scalar>::is_complex_v && vector[1] != 0) {
        throw std::invalid_argument(
            "The distance vector must not have a y-component if the scalar type is real.");
    }

    ion_distance_vector = vector;

    return *this;
}

template <typename Scalar>
SystemAtom<Scalar> &SystemAtom<Scalar>::set_ion_charge(real_t charge) {
    this->hamiltonian_requires_construction = true;
    ion_charge = charge;
    return *this;
}

template <typename Scalar>
SystemAtom<Scalar> &SystemAtom<Scalar>::set_ion_interaction_order(int value) {
    this->hamiltonian_requires_construction = true;
    if (value < 2 || value > 3) {
        throw std::invalid_argument("The order of the Rydberg-ion interaction must be 2 or 3");
    }
    ion_interaction_order = value;
    return *this;
}

template <typename Scalar>
SystemAtom<Scalar> &SystemAtom<Scalar>::set_green_tensor_interpolator(
    const std::shared_ptr<const GreenTensorInterpolator<Scalar>> &green_tensor_interpolator) {
    this->hamiltonian_requires_construction = true;
    this->green_tensor_interpolator = green_tensor_interpolator;
    return *this;
}

template <typename Scalar>
void SystemAtom<Scalar>::construct_hamiltonian() const {
    auto get_operator_matrix = [this](OperatorType type, int q = 0) {
        return this->basis->get_database().get_matrix_elements_in_canonical_basis(
            this->basis, this->basis, type, q);
    };
    // Helper function for constructing matrices of spherical harmonics operators in the
    // canonical basis. The full Hamiltonian is transformed into the actual basis at the end.
    auto get_matrices = [](auto basis, OperatorType type, std::initializer_list<int> m,
                           bool conjugate) {
        std::vector<Eigen::SparseMatrix<Scalar, Eigen::RowMajor>> matrices;
        matrices.reserve(m.size());
        int factor = conjugate ? -1 : 1;
        std::transform(m.begin(), m.end(), std::back_inserter(matrices), [&](int q) {
            return (std::pow(factor, q) *
                    basis->get_database().get_matrix_elements_in_canonical_basis(basis, basis, type,
                                                                                 factor * q))
                .eval();
        });
        return matrices;
    };

    // Construct the unperturbed Hamiltonian in the canonical atomic basis
    this->matrix = utils::get_energies_in_canonical_basis(this->basis);

    this->hamiltonian_is_diagonal = false;
    bool sort_by_quantum_number_f = this->basis->has_quantum_number_f();
    bool sort_by_quantum_number_m = this->basis->has_quantum_number_m();
    bool sort_by_parity = this->basis->has_parity();

    // Estimate the numerical precision so that we can decide which terms to keep
    Eigen::VectorX<real_t> diag = this->matrix.diagonal().real();
    real_t scale = (diag - diag.mean() * Eigen::VectorX<real_t>::Ones(diag.size())).norm();
    real_t numerical_precision = 100 * scale * std::numeric_limits<real_t>::epsilon();

    real_t typical_magnetic_dipole = 1e2;     // ~n^1
    real_t typical_electric_dipole = 1e4;     // ~n^2
    real_t typical_electric_quadrupole = 1e8; // ~n^4

    // Add the interaction with the field of an ion (see
    // https://en.wikipedia.org/wiki/Multipole_expansion#Spherical_form for details)

    Eigen::Map<const Eigen::Vector3<real_t>> vector_map(ion_distance_vector.data(),
                                                        ion_distance_vector.size());
    real_t distance = vector_map.norm();
    SPDLOG_DEBUG("Distance to the ion: {}", distance);

    // Dipole order
    if (std::isfinite(distance) && ion_interaction_order >= 2) {
        // Calculate sqrt(4pi/3) * r * Y_1,q with q = -1, 0, 1 for the first three elements. Take
        // the conjugate of the result and scale it by 1/distance^2.
        Eigen::Vector3<Scalar> vector_dipole_order =
            spherical::get_transformator<Scalar>(1) * vector_map / distance;
        vector_dipole_order = vector_dipole_order.conjugate() / std::pow(distance, 2);

        for (int q = -1; q <= 1; ++q) {
            if (std::abs(vector_dipole_order[q + 1]) * typical_electric_dipole >
                numerical_precision) {
                this->matrix -= ion_charge * vector_dipole_order[q + 1] *
                    get_operator_matrix(OperatorType::ELECTRIC_DIPOLE, q);
                sort_by_quantum_number_f = false;
                sort_by_parity = false;
                sort_by_quantum_number_m &= (q == 0);
            }
        }
    }

    // Quadrupole order (the last entry of vector_quadrupole_order would correspond to
    // ELECTRIC_QUADRUPOLE_ZERO and does not contribute to a *traceless* quadrupole)
    if (std::isfinite(distance) && ion_interaction_order >= 3) {
        // Calculate sqrt(4pi/5) * r^2 * Y_2,q / 3 with q = -2, -1, 0, 1, 2 for the first five
        // elements and (x^2 + y^2 + z^2) / 6 as the last element. Take the conjugate of the
        // result and scale it by 1/distance^3.
        Eigen::Vector<Scalar, 6> vector_quadrupole_order = spherical::get_transformator<Scalar>(2) *
            Eigen::KroneckerProduct(vector_map / distance, vector_map / distance);
        vector_quadrupole_order = 3 * vector_quadrupole_order.conjugate() / std::pow(distance, 3);

        for (int q = -2; q <= 2; ++q) {
            if (std::abs(vector_quadrupole_order[q + 2]) * typical_electric_quadrupole >
                numerical_precision) {
                this->matrix -= ion_charge * vector_quadrupole_order[q + 2] *
                    get_operator_matrix(OperatorType::ELECTRIC_QUADRUPOLE, q);
                sort_by_quantum_number_f = false;
                sort_by_quantum_number_m &= (q == 0);
            }
        }
    }

    // Add external fields (see https://arxiv.org/abs/1612.08053 for details)

    // Transform the electric field to spherical coordinates and take the conjugate
    Eigen::Vector3<Scalar> electric_field_spherical = spherical::get_transformator<Scalar>(1) *
        Eigen::Map<const Eigen::Vector3<real_t>>(electric_field.data(), electric_field.size());
    electric_field_spherical = electric_field_spherical.conjugate();

    // Transform the magnetic field to spherical coordinates and take the conjugate
    Eigen::Vector3<Scalar> magnetic_field_spherical = spherical::get_transformator<Scalar>(1) *
        Eigen::Map<const Eigen::Vector3<real_t>>(magnetic_field.data(), magnetic_field.size());
    magnetic_field_spherical = magnetic_field_spherical.conjugate();

    // Stark effect: - \vec{d} \vec{E} = - d_{1,0} E_{0} + d_{1,1} E_{-} + d_{1,-1} E_{+}
    // = - d_{1,0} E_{0} - d_{1,1} E_{+}^* - d_{1,-1} E_{-}^*
    // with the electric dipole operator: d_{1,q} = - e r sqrt{4 pi / 3} Y_{1,q}(\theta, \phi)
    // where electric_field_spherical=[E_{-}^*, E_{0}, E_{+}^*]
    for (int q = -1; q <= 1; ++q) {
        if (std::abs(electric_field_spherical[q + 1]) * typical_electric_dipole >
            numerical_precision) {
            this->matrix -= electric_field_spherical[q + 1] *
                get_operator_matrix(OperatorType::ELECTRIC_DIPOLE, q);
            sort_by_quantum_number_f = false;
            sort_by_parity = false;
            sort_by_quantum_number_m &= (q == 0);
        }
    }

    // Zeeman effect: - \vec{\mu} \vec{B} = - \mu_{1,0} B_{0} + \mu_{1,1} B_{-} + \mu_{1,-1} B_{+}
    // = - \mu_{1,0} B_{0} - \mu_{1,1} B_{+}^* - \mu_{1,-1} B_{-}^*
    // with the magnetic dipole operator: \vec{\mu} = - \mu_B / \hbar (g_l \vec{l} + g_s \vec{s})
    // where magnetic_field_spherical=[B_{-}^*, B_{0}, B_{+}^*]
    for (int q = -1; q <= 1; ++q) {
        if (std::abs(magnetic_field_spherical[q + 1]) * typical_magnetic_dipole >
            numerical_precision) {
            this->matrix -= magnetic_field_spherical[q + 1] *
                get_operator_matrix(OperatorType::MAGNETIC_DIPOLE, q);
            sort_by_quantum_number_f = false;
            sort_by_quantum_number_m &= (q == 0);
        }
    }

    // Diamagnetism: 1 / (8 m_e) abs(\vec{d} \times \vec{B})^2
    // = (1/12) \left[ B_0^2 (q0 - d_{2,0}) +  B_+ B_- (-2 q0 - d_{2,0})
    // + \sqrt{3} B_0 B_- d_{2,1} + \sqrt{3} B_0 B_+ d_{2,-1}
    // - \sqrt{3/2} B_-^2 d_{2,2} - \sqrt{3/2} B_+^2 d_{2,-2} \right]
    // with the operator: q0 = e^2 r^2 sqrt{4 pi} Y_{0,0} = e^2 r^2
    // and the electric quadrupole operator: d_{2,q} = e^2 r^2 sqrt{4 pi / 5} Y_{2,q}(\theta, \phi)
    // where magnetic_field_spherical=[B_{-}^*, B_{0}, B_{+}^*]
    if (diamagnetism_enabled) {
        if (std::abs(magnetic_field_spherical[1]) * typical_electric_quadrupole >
            numerical_precision) {
            this->matrix += static_cast<real_t>(1 / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[1], 2)) *
                get_operator_matrix(OperatorType::ELECTRIC_QUADRUPOLE_ZERO, 0);
            this->matrix -= static_cast<real_t>(1 / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[1], 2)) *
                get_operator_matrix(OperatorType::ELECTRIC_QUADRUPOLE, 0);
            sort_by_quantum_number_f = false;
        }
        if (std::abs(magnetic_field_spherical[0]) * typical_electric_quadrupole >
                numerical_precision &&
            std::abs(magnetic_field_spherical[2]) * typical_electric_quadrupole >
                numerical_precision) {
            this->matrix -= static_cast<real_t>(2 / 12.) * magnetic_field_spherical[0] *
                magnetic_field_spherical[2] *
                get_operator_matrix(OperatorType::ELECTRIC_QUADRUPOLE_ZERO, 0);
            this->matrix -= static_cast<real_t>(1 / 12.) * magnetic_field_spherical[0] *
                magnetic_field_spherical[2] *
                get_operator_matrix(OperatorType::ELECTRIC_QUADRUPOLE, 0);
            sort_by_quantum_number_f = false;
        }
        if (std::abs(magnetic_field_spherical[1]) * typical_electric_quadrupole >
                numerical_precision &&
            std::abs(magnetic_field_spherical[2]) * typical_electric_quadrupole >
                numerical_precision) {
            this->matrix -= static_cast<real_t>(std::sqrt(3.0) / 12.) *
                magnetic_field_spherical[1] * magnetic_field_spherical[2] *
                get_operator_matrix(OperatorType::ELECTRIC_QUADRUPOLE, 1);
            sort_by_quantum_number_f = false;
            sort_by_quantum_number_m = false;
        }
        if (std::abs(magnetic_field_spherical[1]) * typical_electric_quadrupole >
                numerical_precision &&
            std::abs(magnetic_field_spherical[0]) * typical_electric_quadrupole >
                numerical_precision) {
            this->matrix -= static_cast<real_t>(std::sqrt(3.0) / 12.) *
                magnetic_field_spherical[1] * magnetic_field_spherical[0] *
                get_operator_matrix(OperatorType::ELECTRIC_QUADRUPOLE, -1);
            sort_by_quantum_number_f = false;
            sort_by_quantum_number_m = false;
        }
        if (std::abs(magnetic_field_spherical[2]) * typical_electric_quadrupole >
            numerical_precision) {
            this->matrix -= static_cast<real_t>(std::sqrt(1.5) / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[2], 2)) *
                get_operator_matrix(OperatorType::ELECTRIC_QUADRUPOLE, 2);
            sort_by_quantum_number_f = false;
            sort_by_quantum_number_m = false;
        }
        if (std::abs(magnetic_field_spherical[0]) * typical_electric_quadrupole >
            numerical_precision) {
            this->matrix -= static_cast<real_t>(std::sqrt(1.5) / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[0], 2)) *
                get_operator_matrix(OperatorType::ELECTRIC_QUADRUPOLE, -2);
            sort_by_quantum_number_f = false;
            sort_by_quantum_number_m = false;
        }
    }

    // Add self-interaction via Green tensor
    // H_self = (1/2) * Σ_{ij} D_left[i] * G_{ij} * D_right[j]
    // where D_left uses conjugated convention and
    // D_right uses normal convention.

    // Helper function for adding self interaction
    auto add_interaction = [this, &sort_by_quantum_number_f,
                            &sort_by_quantum_number_m](const auto &entries, const auto &op_left,
                                                       const auto &op_right, int delta) {
        for (const auto &entry : entries) {
            if (std::holds_alternative<
                    typename GreenTensorInterpolator<Scalar>::OmegaDependentEntry>(entry)) {
                throw std::logic_error(
                    "Green tensor with omega dependent entries is currently not supported.");
            }

            const auto &constant_entry =
                std::get<typename GreenTensorInterpolator<Scalar>::ConstantEntry>(entry);
            Eigen::SparseMatrix<Scalar, Eigen::RowMajor> self_interaction =
                (constant_entry.val() / 2. * op_left[constant_entry.row()] *
                 op_right[constant_entry.col()])
                    .eval();
            this->matrix += self_interaction;

            sort_by_quantum_number_f = false;
            if (constant_entry.row() != constant_entry.col() + delta) {
                sort_by_quantum_number_m = false;
            }
        }
    };

    if (green_tensor_interpolator) {
        if (!green_tensor_interpolator->get_spherical_entries(1, 1).empty()) {
            auto dipole_left =
                get_matrices(this->basis, OperatorType::ELECTRIC_DIPOLE, {-1, 0, +1}, true);
            auto dipole_right =
                get_matrices(this->basis, OperatorType::ELECTRIC_DIPOLE, {-1, 0, +1}, false);

            add_interaction(green_tensor_interpolator->get_spherical_entries(1, 1), dipole_left,
                            dipole_right, 0);
        }
    }

    // Transform from the canonical basis into the actual basis
    this->matrix =
        this->basis->get_coefficients().adjoint() * this->matrix * this->basis->get_coefficients();

    // Store which labels can be used to block-diagonalize the Hamiltonian
    this->blockdiagonalizing_labels.clear();
    if (sort_by_quantum_number_f) {
        this->blockdiagonalizing_labels.push_back(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
    }
    if (sort_by_quantum_number_m) {
        this->blockdiagonalizing_labels.push_back(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    }
    if (sort_by_parity) {
        this->blockdiagonalizing_labels.push_back(TransformationType::SORT_BY_PARITY);
    }
}

// Explicit instantiations
template class SystemAtom<double>;
template class SystemAtom<std::complex<double>>;
} // namespace pairinteraction
