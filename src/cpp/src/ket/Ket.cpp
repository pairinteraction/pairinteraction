// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/Ket.hpp"

#include "pairinteraction/ket/QuantumNumberNotAvailableError.hpp"
#include "pairinteraction/utils/hash.hpp"

namespace pairinteraction {
Ket::Ket(double energy, std::unordered_map<std::string, double> quantum_numbers)
    : energy(energy), quantum_numbers(std::move(quantum_numbers)) {}

double Ket::get_energy() const { return energy; }

bool Ket::has_quantum_number(const std::string &name) const {
    return quantum_numbers.contains(name);
}

double Ket::get_quantum_number(const std::string &name) const {
    auto it = quantum_numbers.find(name);
    if (it == quantum_numbers.end()) {
        throw QuantumNumberNotAvailableError(name);
    }
    return it->second;
}

bool Ket::operator==(const Ket &other) const {
    return energy == other.energy && quantum_numbers == other.quantum_numbers;
}

size_t Ket::hash::operator()(const Ket &k) const {
    size_t seed = 0;
    utils::hash_combine(seed, k.energy);
    // The quantum numbers are stored in an unordered map, so we combine the per-entry hashes in an
    // order-independent way (via xor) to obtain a deterministic result.
    size_t quantum_numbers_hash = 0;
    for (const auto &[key, value] : k.quantum_numbers) {
        size_t entry_seed = 0;
        utils::hash_combine(entry_seed, key);
        utils::hash_combine(entry_seed, value);
        quantum_numbers_hash ^= entry_seed;
    }
    utils::hash_combine(seed, quantum_numbers_hash);
    return seed;
}
} // namespace pairinteraction
