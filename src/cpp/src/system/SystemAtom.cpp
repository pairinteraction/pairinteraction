#include "pairinteraction/system/SystemAtom.hpp"

#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/utils/spherical.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <set>

namespace pairinteraction {
template <typename Scalar>
SystemAtom<Scalar>::SystemAtom(std::shared_ptr<const basis_t> basis)
    : System<SystemAtom<Scalar>>(std::move(basis)) {}

template <typename Scalar>
SystemAtom<Scalar> &SystemAtom<Scalar>::set_electric_field(const std::array<real_t, 3> &field) {
    this->hamiltonian_requires_construction = true;
    electric_field_spherical = spherical::convert_to_spherical_basis<Scalar>(field);
    return *this;
}

template <typename Scalar>
SystemAtom<Scalar> &SystemAtom<Scalar>::set_magnetic_field(const std::array<real_t, 3> &field) {
    this->hamiltonian_requires_construction = true;
    magnetic_field_spherical = spherical::convert_to_spherical_basis<Scalar>(field);
    return *this;
}

template <typename Scalar>
SystemAtom<Scalar> &SystemAtom<Scalar>::enable_diamagnetism(bool enable) {
    this->hamiltonian_requires_construction = true;
    diamagnetism_enabled = enable;
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

    // Add external fields (see https://arxiv.org/abs/1612.08053 for details)

    // Stark effect: - \vec{d} \vec{E} = - d_{1,0} E_{0} + d_{1,1} E_{-} + d_{1,-1} E_{+}
    // with the electric dipole operator: d_{1,q} = - e r sqrt{4 pi / 3} Y_{1,q}(\theta, \phi)
    // where electric_field_spherical=[E_{-}, E_{0}, E_{+}]
    if (std::abs(electric_field_spherical[1]) > numerical_precision) {
        *this->hamiltonian -= electric_field_spherical[1] *
            OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_DIPOLE, 0);
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
        sort_by_parity = false;
    }
    if (std::abs(electric_field_spherical[2]) > numerical_precision) {
        *this->hamiltonian += electric_field_spherical[2] *
            OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_DIPOLE, -1);
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
        sort_by_quantum_number_m = false;
        sort_by_parity = false;
    }
    if (std::abs(electric_field_spherical[0]) > numerical_precision) {
        *this->hamiltonian += electric_field_spherical[0] *
            OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_DIPOLE, 1);
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
        sort_by_quantum_number_m = false;
        sort_by_parity = false;
    }

    // Zeeman effect: - \vec{\mu} \vec{B} = - \mu_{1,0} B_{0} + \mu_{1,1} B_{-} + \mu_{1,-1} B_{+}
    // with the magnetic dipole operator: \vec{\mu} = - \mu_B / \hbar (g_l \vec{l} + g_s \vec{s})
    // where magnetic_field_spherical=[B_{-}, B_{0}, B_{+}]
    if (std::abs(magnetic_field_spherical[1]) > numerical_precision) {
        *this->hamiltonian -= magnetic_field_spherical[1] *
            OperatorAtom<Scalar>(basis, OperatorType::MAGNETIC_DIPOLE, 0);
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
    }
    if (std::abs(magnetic_field_spherical[2]) > numerical_precision) {
        *this->hamiltonian += magnetic_field_spherical[2] *
            OperatorAtom<Scalar>(basis, OperatorType::MAGNETIC_DIPOLE, -1);
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
        sort_by_quantum_number_m = false;
    }
    if (std::abs(magnetic_field_spherical[0]) > numerical_precision) {
        *this->hamiltonian += magnetic_field_spherical[0] *
            OperatorAtom<Scalar>(basis, OperatorType::MAGNETIC_DIPOLE, 1);
        this->hamiltonian_is_diagonal = false;
        sort_by_quantum_number_f = false;
        sort_by_quantum_number_m = false;
    }

    // Diamagnetism: 1 / (8 m_e) abs(\vec{d} \times \vec{B})^2
    // = (1/12) \left[ B_0^2 (q0 - d_{2,0}) +  B_+ B_- (-2 q0 - d_{2,0})
    // + B_0 B_- (\sqrt{3} d_{2,1}) + B_0 B_+ (\sqrt{3} d_{2,-1})
    // - B_-^2 d_{2,2} - B_+^2 d_{2,-2} \right]
    // with the operator: q0 = e^2 r^2 sqrt{4 pi} Y_{0,0} = e^2 r^2
    // and the electric quadrupole operator: d_{2,q} = e^2 r^2 sqrt{4 pi / 5} Y_{2,q}(\theta, \phi)
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
        if (std::abs(magnetic_field_spherical[2]) > numerical_precision &&
            std::abs(magnetic_field_spherical[0]) > numerical_precision) {
            *this->hamiltonian -= static_cast<real_t>(2 / 12.) * magnetic_field_spherical[2] *
                magnetic_field_spherical[0] *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE_ZERO, 0);
            *this->hamiltonian -= static_cast<real_t>(1 / 12.) * magnetic_field_spherical[2] *
                magnetic_field_spherical[0] *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, 0);
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
        }
        if (std::abs(magnetic_field_spherical[1]) > numerical_precision &&
            std::abs(magnetic_field_spherical[0]) > numerical_precision) {
            *this->hamiltonian += static_cast<real_t>(std::sqrt(3.0) / 12.) *
                magnetic_field_spherical[1] * magnetic_field_spherical[0] *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, 1);
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
            sort_by_quantum_number_m = false;
        }
        if (std::abs(magnetic_field_spherical[1]) > numerical_precision &&
            std::abs(magnetic_field_spherical[2]) > numerical_precision) {
            *this->hamiltonian += static_cast<real_t>(std::sqrt(3.0) / 12.) *
                magnetic_field_spherical[1] * magnetic_field_spherical[2] *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, -1);
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
            sort_by_quantum_number_m = false;
        }
        if (std::abs(magnetic_field_spherical[0]) > numerical_precision) {
            *this->hamiltonian -= static_cast<real_t>(std::sqrt(1.5) / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[0], 2)) *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, 2);
            this->hamiltonian_is_diagonal = false;
            sort_by_quantum_number_f = false;
            sort_by_quantum_number_m = false;
        }
        if (std::abs(magnetic_field_spherical[2]) > numerical_precision) {
            *this->hamiltonian -= static_cast<real_t>(std::sqrt(1.5) / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[2], 2)) *
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
template class SystemAtom<float>;
template class SystemAtom<double>;
template class SystemAtom<std::complex<float>>;
template class SystemAtom<std::complex<double>>;
} // namespace pairinteraction
