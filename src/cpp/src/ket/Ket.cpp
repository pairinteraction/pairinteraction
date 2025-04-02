// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/Ket.hpp"

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <limits>

namespace pairinteraction {
Ket::Ket(double energy, double f, double m, Parity p)
    : energy(energy), quantum_number_f(f), quantum_number_m(m), parity(p) {}

bool Ket::has_quantum_number_f() const {
    return quantum_number_f != std::numeric_limits<double>::max();
}

bool Ket::has_quantum_number_m() const {
    return quantum_number_m != std::numeric_limits<double>::max();
}

bool Ket::has_parity() const { return parity != Parity::UNKNOWN; }

double Ket::get_energy() const { return energy; }

double Ket::get_quantum_number_f() const { return quantum_number_f; }

double Ket::get_quantum_number_m() const { return quantum_number_m; }

Parity Ket::get_parity() const { return parity; }

bool Ket::operator==(const Ket &other) const {
    return energy == other.energy && quantum_number_f == other.quantum_number_f &&
        quantum_number_m == other.quantum_number_m && parity == other.parity;
}

size_t Ket::hash::operator()(const Ket &k) const {
    size_t seed = 0;
    utils::hash_combine(seed, k.energy);
    utils::hash_combine(seed, k.quantum_number_f);
    utils::hash_combine(seed, k.quantum_number_m);
    utils::hash_combine(seed, static_cast<int>(k.parity));
    return seed;
}
} // namespace pairinteraction
