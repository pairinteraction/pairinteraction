// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/KetAtom.hpp"

#include "pairinteraction/utils/hash.hpp"

#include <array>
#include <cctype>
#include <cmath>
#include <fmt/core.h>
#include <fmt/format.h>
#include <string>
#include <string_view>
#include <vector>

namespace pairinteraction {

constexpr std::array<std::string_view, 6> quantum_number_l_labels = {"S", "P", "D", "F", "G", "H"};

KetAtom::KetAtom(Private /*unused*/, double energy, std::string species,
                 std::unordered_map<std::string, double> quantum_numbers, Database &database,
                 size_t id_in_database)
    : Ket(energy), species(std::move(species)), quantum_numbers(std::move(quantum_numbers)),
      database(database), id_in_database(id_in_database) {}

Database &KetAtom::get_database() const { return database; }

size_t KetAtom::get_id_in_database() const { return id_in_database; }

double KetAtom::get_quantum_number(const std::string &name) const {
    if (auto it = quantum_numbers.find("exp_" + name); it != quantum_numbers.end()) {
        return it->second;
    }
    return quantum_numbers.at(name);
}

double KetAtom::get_quantum_number_std(const std::string &name) const {
    if (auto it = quantum_numbers.find("std_" + name); it != quantum_numbers.end()) {
        return it->second;
    }
    return 0;
}

std::string KetAtom::get_label() const {
    double s = get_quantum_number("s");
    double l = get_quantum_number("l");
    double f = get_quantum_number("f");
    double m = get_quantum_number("m");
    bool is_calculated_with_mqdt = get_quantum_number("is_calculated_with_mqdt") != 0;

    size_t pos = species.find('_');
    std::string label = (pos != std::string::npos) ? species.substr(0, pos) : species;
    label[0] = static_cast<char>(std::toupper(label[0]));

    if (!is_calculated_with_mqdt) {
        if (s == 0) {
            label += "_singlet";
        } else if (s == 1) {
            label += "_triplet";
        } else if (s != 0.5) {
            throw std::runtime_error(
                "Invalid value for quantum number s in the single-channel description.");
        }
    }

    label += ":";

    if (is_calculated_with_mqdt) {
        label += fmt::format("S={:.1f},nu={:.1f},L={:.1f},", s, get_quantum_number("nu"), l);
        label += get_quantum_number("is_j_total_momentum") != 0 ? "J=" : "F=";
    } else {
        label += fmt::format("{:.0f},", get_quantum_number("n"));
        if (l == std::rint(l) && l < quantum_number_l_labels.size()) {
            label += quantum_number_l_labels.at(static_cast<size_t>(l));
        } else {
            label += fmt::format("{:.0f}", l);
        }
        label += "_";
    }

    if (f == std::rint(f)) {
        label += fmt::format("{:.0f}", f);
    } else if (2 * f == std::rint(2 * f)) {
        label += fmt::format("{:.0f}/2", 2 * f);
    } else {
        std::abort(); // can't happen because the total momentum is validated to be an integer
                      // or half-integer
    }

    if (m == std::rint(m)) {
        label += fmt::format(",{:.0f}", m);
    } else if (2 * m == std::rint(2 * m)) {
        label += fmt::format(",{:.0f}/2", 2 * m);
    } else {
        std::abort(); // can't happen because the quantum number m is validated to be an integer
                      // or half-integer
    }

    return label;
}

std::shared_ptr<KetAtom>
KetAtom::get_ket_for_different_quantum_number_m(double new_quantum_number_m) const {
    auto ket = *this;
    ket.quantum_numbers.at("m") = new_quantum_number_m;
    return std::make_shared<KetAtom>(ket);
}

const std::string &KetAtom::get_species() const { return species; }

bool KetAtom::operator==(const KetAtom &other) const {
    return Ket::operator==(other) && species == other.species &&
        quantum_numbers == other.quantum_numbers;
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
    utils::hash_combine(seed, quantum_numbers_hash);
    return seed;
}

} // namespace pairinteraction
