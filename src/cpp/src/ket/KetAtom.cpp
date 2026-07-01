// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/KetAtom.hpp"

#include "pairinteraction/utils/hash.hpp"

#include <string>
#include <vector>

namespace pairinteraction {

KetAtom::KetAtom(Private /*unused*/, double energy, std::string species,
                 std::unordered_map<std::string, double> quantum_numbers,
                 std::unordered_map<std::string, double> quantum_numbers_std, Database &database,
                 size_t id_in_database)
    : Ket(energy), species(std::move(species)), quantum_numbers(std::move(quantum_numbers)),
      quantum_numbers_std(std::move(quantum_numbers_std)), database(database),
      id_in_database(id_in_database) {}

Database &KetAtom::get_database() const { return database; }

size_t KetAtom::get_id_in_database() const { return id_in_database; }

double KetAtom::get_quantum_number(const std::string &name) const {
    return quantum_numbers.at(name);
}

double KetAtom::get_quantum_number_std(const std::string &name) const {
    auto it = quantum_numbers_std.find(name);
    return it != quantum_numbers_std.end() ? it->second : 0;
}

const std::string &KetAtom::get_species() const { return species; }

bool KetAtom::operator==(const KetAtom &other) const {
    return Ket::operator==(other) && species == other.species &&
        quantum_numbers == other.quantum_numbers &&
        quantum_numbers_std == other.quantum_numbers_std;
}

bool KetAtom::operator!=(const KetAtom &other) const { return !(*this == other); }

size_t KetAtom::hash::operator()(const KetAtom &k) const {
    size_t seed = typename Ket::hash()(k);
    utils::hash_combine(seed, k.species);
    // The quantum numbers are stored in an unordered map, so we combine the per-entry hashes in an
    // order-independent way (via xor) to obtain a deterministic result.
    size_t quantum_numbers_hash = 0;
    for (const auto &[key, value] : k.quantum_numbers) {
        size_t entry_seed = 0;
        utils::hash_combine(entry_seed, key);
        utils::hash_combine(entry_seed, value);
        quantum_numbers_hash ^= entry_seed;
    }
    for (const auto &[key, value] : k.quantum_numbers_std) {
        size_t entry_seed = 0;
        utils::hash_combine(entry_seed, key);
        utils::hash_combine(entry_seed, value);
        quantum_numbers_hash ^= entry_seed;
    }
    utils::hash_combine(seed, quantum_numbers_hash);
    return seed;
}

} // namespace pairinteraction
