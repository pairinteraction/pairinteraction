#include "pairinteraction/system/SystemAtom.hpp"

#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/utils/spherical.hpp"

#include <set>

namespace pairinteraction {
template <typename Scalar>
SystemAtom<Scalar>::SystemAtom(std::shared_ptr<const basis_t> basis)
    : System<SystemAtom<Scalar>>(basis) {}

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
    this->hamiltonian = std::make_unique<OperatorAtom<Scalar>>(basis, OperatorType::ENERGY);
    this->blockdiagonalizing_labels = {TransformationType::SORT_BY_QUANTUM_NUMBER_F,
                                       TransformationType::SORT_BY_QUANTUM_NUMBER_M,
                                       TransformationType::SORT_BY_PARITY};

    real_t precision = 10 * std::numeric_limits<real_t>::epsilon();

    // Stark effect
    if (std::abs(electric_field_spherical[0]) > precision) {
        *this->hamiltonian -= electric_field_spherical[0] *
            OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_DIPOLE, 0);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_PARITY);
    }
    if (std::abs(electric_field_spherical[1]) > precision) {
        *this->hamiltonian += electric_field_spherical[1] *
            OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_DIPOLE, -1);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_PARITY);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    }
    if (std::abs(electric_field_spherical[2]) > precision) {
        *this->hamiltonian += electric_field_spherical[2] *
            OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_DIPOLE, 1);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_PARITY);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    }

    // Zeeman effect
    if (std::abs(magnetic_field_spherical[0]) > precision) {
        *this->hamiltonian -= magnetic_field_spherical[0] *
            OperatorAtom<Scalar>(basis, OperatorType::MAGNETIC_DIPOLE, 0);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
    }
    if (std::abs(magnetic_field_spherical[1]) > precision) {
        *this->hamiltonian += magnetic_field_spherical[1] *
            OperatorAtom<Scalar>(basis, OperatorType::MAGNETIC_DIPOLE, -1);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    }
    if (std::abs(magnetic_field_spherical[2]) > precision) {
        *this->hamiltonian += magnetic_field_spherical[2] *
            OperatorAtom<Scalar>(basis, OperatorType::MAGNETIC_DIPOLE, 1);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
        this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    }

    // Diamagnetism
    if (diamagnetism_enabled) {
        if (std::abs(magnetic_field_spherical[0]) > precision) {
            *this->hamiltonian += static_cast<real_t>(1 / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[0], 2)) *
                OperatorAtom<Scalar>(basis, OperatorType::DIAMAGNETIC, 0);
            *this->hamiltonian -= static_cast<real_t>(1 / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[0], 2)) *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, 0);
            this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
        }
        if (std::abs(magnetic_field_spherical[1]) > precision &&
            std::abs(magnetic_field_spherical[2]) > precision) {
            *this->hamiltonian -= static_cast<real_t>(2 / 12.) * magnetic_field_spherical[1] *
                magnetic_field_spherical[2] *
                OperatorAtom<Scalar>(basis, OperatorType::DIAMAGNETIC, 0);
            *this->hamiltonian -= static_cast<real_t>(1 / 12.) * magnetic_field_spherical[1] *
                magnetic_field_spherical[2] *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, 0);
            this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
            this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
        }
        if (std::abs(magnetic_field_spherical[0]) > precision &&
            std::abs(magnetic_field_spherical[2]) > precision) {
            *this->hamiltonian += static_cast<real_t>(std::sqrt(3.0) / 12.) *
                magnetic_field_spherical[0] * magnetic_field_spherical[2] *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, 1);
            this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
            this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
        }
        if (std::abs(magnetic_field_spherical[0]) > precision &&
            std::abs(magnetic_field_spherical[1]) > precision) {
            *this->hamiltonian += static_cast<real_t>(std::sqrt(3.0) / 12.) *
                magnetic_field_spherical[0] * magnetic_field_spherical[1] *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, -1);
            this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
            this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
        }
        if (std::abs(magnetic_field_spherical[2]) > precision) {
            *this->hamiltonian -= static_cast<real_t>(std::sqrt(1.5) / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[2], 2)) *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, 2);
            this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
            this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
        }
        if (std::abs(magnetic_field_spherical[1]) > precision) {
            *this->hamiltonian -= static_cast<real_t>(std::sqrt(1.5) / 12.) *
                static_cast<Scalar>(std::pow(magnetic_field_spherical[1], 2)) *
                OperatorAtom<Scalar>(basis, OperatorType::ELECTRIC_QUADRUPOLE, -2);
            this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_F);
            this->blockdiagonalizing_labels.erase(TransformationType::SORT_BY_QUANTUM_NUMBER_M);
        }
    }
}

// Explicit instantiations
template class SystemAtom<float>;
template class SystemAtom<double>;
template class SystemAtom<std::complex<float>>;
template class SystemAtom<std::complex<double>>;
} // namespace pairinteraction
