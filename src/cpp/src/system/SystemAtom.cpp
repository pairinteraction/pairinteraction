// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/SystemAtom.hpp"

#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/spherical.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <limits>
#include <memory>
#include <set>
#include <spdlog/spdlog.h>

namespace pairinteraction {
template <typename Scalar>
SystemAtom<Scalar>::SystemAtom(std::shared_ptr<const basis_t> basis)
    : System<SystemAtom<Scalar>>(std::move(basis)) {}

template <typename Scalar>
SystemAtom<Scalar> &SystemAtom<Scalar>::set_electric_field(const std::array<real_t, 3> &field) {
    this->hamiltonian_requires_construction = true;

    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();
    if (!traits::NumTraits<Scalar>::is_complex_v && std::abs(field[1]) > numerical_precision) {
        throw std::invalid_argument(
            "The field must not have a y-component if the scalar type is real.");
    }

    Eigen::Map<const Eigen::Vector3<real_t>> field_map(field.data(), field.size());
    Eigen::Map<Eigen::Vector3<Scalar>> electric_field_spherical_map(
        electric_field_spherical.data(), electric_field_spherical.size());
    electric_field_spherical_map = spherical::get_transformator<Scalar>(1) * field_map;
    electric_field_spherical_map = electric_field_spherical_map.conjugate();

    return *this;
}

template <typename Scalar>
SystemAtom<Scalar> &SystemAtom<Scalar>::set_magnetic_field(const std::array<real_t, 3> &field) {
    this->hamiltonian_requires_construction = true;

    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();
    if (!traits::NumTraits<Scalar>::is_complex_v && std::abs(field[1]) > numerical_precision) {
        throw std::invalid_argument(
            "The field must not have a y-component if the scalar type is real.");
    }

    Eigen::Map<const Eigen::Vector3<real_t>> field_map(field.data(), field.size());
    Eigen::Map<Eigen::Vector3<Scalar>> magnetic_field_spherical_map(
        magnetic_field_spherical.data(), magnetic_field_spherical.size());
    magnetic_field_spherical_map = spherical::get_transformator<Scalar>(1) * field_map;
    magnetic_field_spherical_map = magnetic_field_spherical_map.conjugate();

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
SystemAtom<Scalar>::set_distance_vector_to_ion(const std::array<real_t, 3> &vector) {
    this->hamiltonian_requires_construction = true;

    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();
    if (!traits::NumTraits<Scalar>::is_complex_v && std::abs(vector[1]) > numerical_precision) {
        throw std::invalid_argument(
            "The distance vector must not have a y-component if the scalar type is real.");
    }

    distance_vector_to_ion = vector;

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
    if (order < 2 || order > 3) {
        throw std::invalid_argument("The order of the Rydberg-ion interaction must be 2 or 3");
    }
    order = value;
    return *this;
}

template <typename Scalar>
void SystemAtom<Scalar>::construct_hamiltonian() const {
    auto basis = this->hamiltonian->get_basis();

    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

    // Construct the unperturbed Hamiltonian
    this->hamiltonian = std::make_unique<OperatorAtom<Scalar>>(basis, OperatorType::ENERGY);
    this->hamiltonian_is_diagonal = true;
    bool sort_by_quantum_number_f = true;
    bool sort_by_quantum_number_m = true;
    bool sort_by_parity = true;

    // Add the interaction with the field of an ion (see
    // https://en.wikipedia.org/wiki/Multipole_expansion#Spherical_form for details)

    Eigen::Map<const Eigen::Vector3<real_t>> vector_map(distance_vector_to_ion.data(),
                                                        distance_vector_to_ion.size());
    real_t distance = vector_map.norm();
    SPDLOG_DEBUG("Distance to the ion: {}", distance);

    // Dipole order
    if (distance != std::numeric_limits<real_t>::infinity() && order >= 2) {
        // Calculate sqrt(4pi/3) * r * Y_1,q with q = -1, 0, 1 for the first three elements and take
        // the conjugate of the result
        Eigen::Vector3<Scalar> ion_first_order =
            spherical::get_transformator<Scalar>(1) * vector_map / distance;
        ion_first_order = ion_first_order.conjugate() / std::pow(distance, 2);

        for (int q = -1; q <= 1; ++q) {
            if (std::abs(ion_first_order[q + 1]) > numerical_precision) {
                *this->hamiltonian -= ion_charge * ion_first_order[q + 1] *
                    OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_DIPOLE, q);
                this->hamiltonian_is_diagonal = false;
                sort_by_quantum_number_f = false;
                sort_by_parity = false;
                sort_by_quantum_number_m &= (q == 0);
            }
        }
    }

    // Quadrupole order (the last entry of ion_second_order would correspond to
    // ELECTRIC_QUADRUPOLE_ZERO and does not contribute to a *traceless* quadrupole)
    if (distance != std::numeric_limits<real_t>::infinity() && order >= 3) {
        // Calculate sqrt(4pi/5) * r^2 * Y_2,q / 3 with q = -2, -1, 0, 1, 2 for the first five
        // elements and (x^2 + y^2 + z^2) / 6 as the last element, and take the conjugate of the
        // result
        Eigen::Vector<Scalar, 6> ion_second_order = spherical::get_transformator<Scalar>(2) *
            Eigen::KroneckerProduct(vector_map / distance, vector_map / distance);
        ion_second_order = 3 * ion_second_order.conjugate() / std::pow(distance, 3);

        for (int q = -2; q <= 2; ++q) {
            if (std::abs(ion_second_order[q + 2]) > numerical_precision) {
                *this->hamiltonian -= ion_charge * ion_second_order[q + 2] *
                    OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, q);
                this->hamiltonian_is_diagonal = false;
                sort_by_quantum_number_f = false;
                sort_by_quantum_number_m &= (q == 0);
            }
        }
    }

    // Stark effect: - \vec{d} \vec{E} = - d_{1,0} E_{0} + d_{1,1} E_{-} + d_{1,-1} E_{+}
    // = - d_{1,0} E_{0} - d_{1,1} E_{+}^* - d_{1,-1} E_{-}^*
    // with the electric dipole operator: d_{1,q} = - e r sqrt{4 pi / 3} Y_{1,q}(\theta, \phi)
    // where electric_field_spherical=[E_{-}^*, E_{0}, E_{+}^*]

    for (int q = -1; q <= 1; ++q) {
        if (std::abs(electric_field_spherical.at(q + 1)) > numerical_precision) {
            *this->hamiltonian -= electric_field_spherical.at(q + 1) *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_DIPOLE, q);
            this->hamiltonian_is_diagonal = false;
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
        if (std::abs(magnetic_field_spherical.at(q + 1)) > numerical_precision) {
            *this->hamiltonian -= magnetic_field_spherical.at(q + 1) *
                OperatorAtom<Scalar>(basis, OperatorType::MAGNETIC_DIPOLE, q);
            this->hamiltonian_is_diagonal = false;
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
        if (std::abs(magnetic_field_spherical[1]) > numerical_precision) {
            *this->hamiltonian += static_cast<real_t>(1 / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[1], 2)) *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE_ZERO, 0);
            *this->hamiltonian -= static_cast<real_t>(1 / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[1], 2)) *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, 0);
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
        }
        if (std::abs(magnetic_field_spherical[0]) > numerical_precision &&
            std::abs(magnetic_field_spherical[2]) > numerical_precision) {
            *this->hamiltonian -= static_cast<real_t>(2 / 12.) * magnetic_field_spherical[0] *
                magnetic_field_spherical[2] *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE_ZERO, 0);
            *this->hamiltonian -= static_cast<real_t>(1 / 12.) * magnetic_field_spherical[0] *
                magnetic_field_spherical[2] *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, 0);
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
        }
        if (std::abs(magnetic_field_spherical[1]) > numerical_precision &&
            std::abs(magnetic_field_spherical[2]) > numerical_precision) {
            *this->hamiltonian -= static_cast<real_t>(std::sqrt(3.0) / 12.) *
                magnetic_field_spherical[1] * magnetic_field_spherical[2] *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, 1);
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
            sort_by_quantum_number_m = false;
        }
        if (std::abs(magnetic_field_spherical[1]) > numerical_precision &&
            std::abs(magnetic_field_spherical[0]) > numerical_precision) {
            *this->hamiltonian -= static_cast<real_t>(std::sqrt(3.0) / 12.) *
                magnetic_field_spherical[1] * magnetic_field_spherical[0] *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, -1);
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
            sort_by_quantum_number_m = false;
        }
        if (std::abs(magnetic_field_spherical[0]) > numerical_precision) {
            *this->hamiltonian -= static_cast<real_t>(std::sqrt(1.5) / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[0], 2)) *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, -2);
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
            sort_by_quantum_number_m = false;
        }
    }

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
