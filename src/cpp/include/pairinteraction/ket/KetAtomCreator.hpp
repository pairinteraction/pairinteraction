// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/enums/Parity.hpp"

#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>

namespace pairinteraction {
class Database;

class KetAtom;

/**
 * @class KetAtomCreator
 *
 * @brief Builder class for creating KetAtom objects.
 */
class KetAtomCreator {
public:
    KetAtomCreator() = default;
    KetAtomCreator(std::string species, int n, double l, double j, double m);
    KetAtomCreator &set_species(const std::string &value);
    KetAtomCreator &set_energy(double value);
    KetAtomCreator &set_parity(Parity value);
    // Set the quantum number with the given logical name (e.g. "f", "m", "n", "nu", "l", ...).
    KetAtomCreator &set_quantum_number(const std::string &name, double value);
    std::shared_ptr<const KetAtom> create(Database &database) const;

private:
    std::optional<std::string> species;
    Parity parity{Parity::UNKNOWN};
    std::optional<double> energy;
    std::unordered_map<std::string, double> quantum_numbers;
};

} // namespace pairinteraction
