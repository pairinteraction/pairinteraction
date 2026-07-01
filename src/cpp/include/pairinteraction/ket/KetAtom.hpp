// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/ket/Ket.hpp"

#include <string>
#include <type_traits>
#include <unordered_map>

namespace pairinteraction {
class Database;

/**
 * @class KetAtom
 *
 * @brief Class for representing atomic kets.
 */
class KetAtom : public Ket {
    friend class Database;
    struct Private {};

public:
    KetAtom(Private /*unused*/, double energy, std::string species,
            std::unordered_map<std::string, double> quantum_numbers,
            std::unordered_map<std::string, double> quantum_numbers_std, Database &database,
            size_t id_in_database);

    Database &get_database() const;
    size_t get_id_in_database() const;
    const std::string &get_species() const;
    double get_quantum_number(const std::string &name) const;
    double get_quantum_number_std(const std::string &name) const;

    bool operator==(const KetAtom &other) const;
    bool operator!=(const KetAtom &other) const;

    struct hash {
        std::size_t operator()(const KetAtom &k) const;
    };

private:
    std::string species;
    std::unordered_map<std::string, double> quantum_numbers;
    std::unordered_map<std::string, double> quantum_numbers_std;
    Database &database;
    size_t id_in_database;
};

} // namespace pairinteraction
