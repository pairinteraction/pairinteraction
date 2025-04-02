// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/enums/Parity.hpp"

#include <memory>
#include <optional>
#include <string>
#include <type_traits>

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
    KetAtomCreator &set_quantum_number_f(double value);
    KetAtomCreator &set_quantum_number_m(double value);
    KetAtomCreator &set_parity(Parity value);
    KetAtomCreator &set_quantum_number_n(int value);
    KetAtomCreator &set_quantum_number_nu(double value);
    KetAtomCreator &set_quantum_number_nui(double value);
    KetAtomCreator &set_quantum_number_l(double value);
    KetAtomCreator &set_quantum_number_s(double value);
    KetAtomCreator &set_quantum_number_j(double value);
    KetAtomCreator &set_quantum_number_l_ryd(double value);
    KetAtomCreator &set_quantum_number_j_ryd(double value);
    std::shared_ptr<const KetAtom> create(Database &database) const;

private:
    std::optional<std::string> species;
    Parity parity{Parity::UNKNOWN};
    std::optional<double> energy;
    std::optional<double> quantum_number_f;
    std::optional<double> quantum_number_m;
    std::optional<int> quantum_number_n;
    std::optional<double> quantum_number_nu;
    std::optional<double> quantum_number_nui;
    std::optional<double> quantum_number_l;
    std::optional<double> quantum_number_s;
    std::optional<double> quantum_number_j;
    std::optional<double> quantum_number_l_ryd;
    std::optional<double> quantum_number_j_ryd;
};

} // namespace pairinteraction
