// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/KetAtomCreator.hpp"

#include "pairinteraction/database/AtomDescriptionByParameters.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/Parity.hpp"

#include <cmath>

namespace pairinteraction {
KetAtomCreator::KetAtomCreator(std::string species, int n, double l, double j, double m)
    : species(std::move(species)), quantum_number_f(j), quantum_number_m(m), quantum_number_n(n),
      quantum_number_l(l), quantum_number_s(0.5), quantum_number_j(j) {}

KetAtomCreator &KetAtomCreator::set_species(const std::string &value) {
    species.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_energy(double value) {
    energy.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_quantum_number_f(double value) {
    if (2 * value != std::rint(2 * value)) {
        throw std::invalid_argument("Quantum number f must be an integer or half-integer.");
    }
    quantum_number_f.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_quantum_number_m(double value) {
    if (2 * value != std::rint(2 * value)) {
        throw std::invalid_argument("Quantum number m must be an integer or half-integer.");
    }
    quantum_number_m.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_parity(Parity value) {
    parity = value;
    return *this;
}

KetAtomCreator &KetAtomCreator::set_quantum_number_n(int value) {
    quantum_number_n.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_quantum_number_nu(double value) {
    quantum_number_nu.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_quantum_number_nui(double value) {
    quantum_number_nui.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_quantum_number_l(double value) {
    quantum_number_l.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_quantum_number_s(double value) {
    quantum_number_s.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_quantum_number_j(double value) {
    quantum_number_j.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_quantum_number_l_ryd(double value) {
    quantum_number_l_ryd.emplace(value);
    return *this;
}

KetAtomCreator &KetAtomCreator::set_quantum_number_j_ryd(double value) {
    quantum_number_j_ryd.emplace(value);
    return *this;
}

std::shared_ptr<const KetAtom> KetAtomCreator::create(Database &database) const {

    if (!species.has_value()) {
        throw std::runtime_error("Species not set.");
    }

    AtomDescriptionByParameters description{parity,
                                            energy,
                                            quantum_number_f,
                                            quantum_number_m,
                                            quantum_number_n,
                                            quantum_number_nu,
                                            quantum_number_nui,
                                            quantum_number_l,
                                            quantum_number_s,
                                            quantum_number_j,
                                            quantum_number_l_ryd,
                                            quantum_number_j_ryd};

    return database.get_ket(species.value(), description);
}
} // namespace pairinteraction
