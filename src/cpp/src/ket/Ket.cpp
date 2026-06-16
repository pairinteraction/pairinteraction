// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/Ket.hpp"

#include "pairinteraction/utils/hash.hpp"

namespace pairinteraction {
Ket::Ket(double energy) : energy(energy) {}

double Ket::get_energy() const { return energy; }

bool Ket::operator==(const Ket &other) const { return energy == other.energy; }

size_t Ket::hash::operator()(const Ket &k) const {
    size_t seed = 0;
    utils::hash_combine(seed, k.energy);
    return seed;
}
} // namespace pairinteraction
