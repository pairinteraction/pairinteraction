// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisAtomCreator.hpp"

#include "pairinteraction/database/AtomDescriptionByRanges.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/KetAtom.hpp"

namespace pairinteraction {
template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::set_species(const std::string &value) {
    species.emplace(value);
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &BasisAtomCreator<Scalar>::restrict_energy(real_t min, real_t max) {
    range_energy = {min, max};
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &
BasisAtomCreator<Scalar>::restrict_quantum_number(const std::string &name, real_t min, real_t max) {
    if ((name == "f" || name == "m") &&
        (2 * min != std::rint(2 * min) || 2 * max != std::rint(2 * max))) {
        throw std::invalid_argument("Quantum number " + name +
                                    " must be an integer or half-integer.");
    }
    if (name == "parity" && ((min != 1 && min != -1) || (max != 1 && max != -1))) {
        throw std::invalid_argument("The parity must be +1 or -1.");
    }
    quantum_number_ranges[name] = {min, max};
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &
BasisAtomCreator<Scalar>::set_quantum_number_standard_deviation_factor(real_t value) {
    if (value < 0) {
        throw std::invalid_argument(
            "The quantum number standard deviation factor must be non-negative.");
    }
    quantum_number_standard_deviation_factor = value;
    return *this;
}

template <typename Scalar>
BasisAtomCreator<Scalar> &
BasisAtomCreator<Scalar>::add_ket(const std::shared_ptr<const ket_t> &ket) {
    if (additional_ket_species.has_value() &&
        additional_ket_species.value() != ket->get_species()) {
        throw std::invalid_argument("Species mismatch.");
    }
    additional_ket_species.emplace(ket->get_species());
    additional_ket_ids.push_back(ket->get_id_in_database());
    return *this;
}

template <typename Scalar>
std::shared_ptr<const BasisAtom<Scalar>>
BasisAtomCreator<Scalar>::create(Database &database) const {

    if (species.has_value() && additional_ket_species.has_value() &&
        species.value() != additional_ket_species.value()) {
        throw std::invalid_argument("Species mismatch.");
    }

    std::string extracted_species;
    if (species.has_value()) {
        extracted_species = species.value();
    } else if (additional_ket_species.has_value()) {
        extracted_species = additional_ket_species.value();
    } else {
        throw std::runtime_error("Species not set.");
    }

    AtomDescriptionByRanges description{range_energy, quantum_number_ranges,
                                        quantum_number_standard_deviation_factor};

    return database.get_basis<Scalar>(extracted_species, description, additional_ket_ids);
}

// Explicit instantiations
template class BasisAtomCreator<double>;
template class BasisAtomCreator<std::complex<double>>;
} // namespace pairinteraction
