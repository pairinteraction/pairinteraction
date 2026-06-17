// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/KetAtomCreator.hpp"

#include "pairinteraction/database/AtomDescriptionByParameters.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/Parity.hpp"

#include <cmath>

namespace pairinteraction {
KetAtomCreator::KetAtomCreator(std::string species, int n, double l, double j, double m)
    : species(std::move(species)) {
    set_quantum_number("n", n);
    set_quantum_number("l", l);
    set_quantum_number("s", 0.5);
    set_quantum_number("j", j);
    set_quantum_number("f", j);
    set_quantum_number("m", m);
}

KetAtomCreator &KetAtomCreator::set_species(const std::string &value) {
    species.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_energy(double value) {
    energy.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_parity(Parity value) {
    parity = value;
    return *this;
}

KetAtomCreator &KetAtomCreator::set_quantum_number(const std::string &name, double value) {
    if ((name == "f" || name == "m") && 2 * value != std::rint(2 * value)) {
        throw std::invalid_argument("Quantum number " + name +
                                    " must be an integer or half-integer.");
    }
    quantum_numbers[name] = value;
    return *this;
}

std::shared_ptr<const KetAtom> KetAtomCreator::create(Database &database) const {

    if (!species.has_value()) {
        throw std::runtime_error("Species not set.");
    }

    AtomDescriptionByParameters description{parity, energy, quantum_numbers};

    return database.get_ket(species.value(), description);
}
} // namespace pairinteraction
