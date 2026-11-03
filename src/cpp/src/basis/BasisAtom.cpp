// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisAtom.hpp"

#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/KetAtom.hpp"

#include <cassert>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
// NOLINTNEXTLINE(cert-err58-cpp)
const std::map<SorterType, std::string> BasisAtom<Scalar>::sorter_type_to_quantum_number_name{
    {SorterType::QUANTUM_NUMBER_F, "f"},
    {SorterType::QUANTUM_NUMBER_M, "m"},
    {SorterType::PARITY, "parity"},
};

template <typename Scalar>
BasisAtom<Scalar>::BasisAtom(Private /*unused*/, std::shared_ptr<const ketvec_t> kets,
                             std::shared_ptr<const map_size_t> ket_id_to_ket_index,
                             std::string &&canonical_basis_id, Database &database,
                             Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &&coefficients,
                             energies_t &&energy_of_states,
                             quantum_numbers_t &&quantum_numbers_of_states)
    : Basis<BasisAtom<Scalar>>(std::move(kets), std::move(coefficients),
                               std::move(energy_of_states), std::move(quantum_numbers_of_states)),
      canonical_basis_id(std::move(canonical_basis_id)), database(database),
      ket_id_to_ket_index(std::move(ket_id_to_ket_index)) {}

template <typename Scalar>
std::shared_ptr<const BasisAtom<Scalar>>
BasisAtom<Scalar>::get_modified_copy(Eigen::SparseMatrix<Scalar, Eigen::RowMajor> coefficients,
                                     energies_t energy_of_states,
                                     quantum_numbers_t quantum_numbers_of_states) const {
    return std::make_shared<const Type>(
        Private(), this->kets, ket_id_to_ket_index, std::string(canonical_basis_id), database,
        std::move(coefficients), std::move(energy_of_states), std::move(quantum_numbers_of_states));
}

template <typename Scalar>
Database &BasisAtom<Scalar>::get_database() const {
    return database;
}

template <typename Scalar>
const std::string &BasisAtom<Scalar>::get_species() const {
    return this->get_kets()[0]->get_species();
}

template <typename Scalar>
int BasisAtom<Scalar>::get_ket_index_from_id(size_t ket_id) const {
    auto it = ket_id_to_ket_index->find(ket_id);
    if (it == ket_id_to_ket_index->end()) {
        return -1;
    }
    return static_cast<int>(it->second);
}

template <typename Scalar>
const std::string &BasisAtom<Scalar>::get_canonical_basis_id() const {
    return canonical_basis_id;
}

template <typename Scalar>
std::shared_ptr<const typename BasisAtom<Scalar>::Type>
BasisAtom<Scalar>::merge(const std::shared_ptr<const Type> &other) const {
    if (&database != &other->database) {
        throw std::invalid_argument("Cannot merge atomic bases from different Database instances.");
    }
    if (get_species() != other->get_species()) {
        throw std::invalid_argument("Cannot merge atomic bases with different species.");
    }

    // The database creates the merged kets, sorted by their id in the database
    BasisAtomCreator<Scalar> creator(database);
    std::unordered_set<size_t> ket_ids;
    ket_ids.reserve(this->get_number_of_kets() + other->get_number_of_kets());
    for (const auto &basis : {this, other.get()}) {
        for (const auto &ket : basis->get_kets()) {
            if (ket_ids.insert(ket->get_id_in_database()).second) {
                creator.add_ket(ket);
            }
        }
    }
    auto merged = creator.create();

    auto get_ket_indices = [&](const Type &basis) {
        std::vector<size_t> ket_indices;
        ket_indices.reserve(basis.get_number_of_kets());
        for (const auto &ket : basis.get_kets()) {
            int ket_index = merged->get_ket_index_from_id(ket->get_id_in_database());
            if (ket_index < 0) {
                throw std::runtime_error("The merged basis is missing a ket of the input bases.");
            }
            ket_indices.push_back(static_cast<size_t>(ket_index));
        }
        return ket_indices;
    };
    auto [coefficients, energy_of_states, quantum_numbers_of_states] = this->get_merged_states(
        *other, get_ket_indices(*this), get_ket_indices(*other), merged->get_number_of_kets());
    return merged->get_modified_copy(std::move(coefficients), std::move(energy_of_states),
                                     std::move(quantum_numbers_of_states));
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
BasisAtom<Scalar>::get_matrix_elements(const std::shared_ptr<const Type> &other, OperatorType type,
                                       int q) const {
    auto matrix_elements = this->get_database().get_matrix_elements_in_canonical_basis(
        this->shared_from_this(), other, type, q);
    matrix_elements =
        other->get_coefficients().adjoint() * matrix_elements * this->get_coefficients();

    assert(static_cast<size_t>(matrix_elements.rows()) == other->get_number_of_states());
    assert(static_cast<size_t>(matrix_elements.cols()) == this->get_number_of_states());

    return matrix_elements;
}

// Explicit instantiations
template class BasisAtom<double>;
template class BasisAtom<std::complex<double>>;
} // namespace pairinteraction
