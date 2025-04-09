// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/KetAtom.hpp"

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <array>
#include <cctype>
#include <cmath>
#include <fmt/core.h>
#include <string>
#include <string_view>
#include <vector>

namespace pairinteraction {

constexpr std::array<std::string_view, 6> quantum_number_l_labels = {"S", "P", "D", "F", "G", "H"};

KetAtom::KetAtom(Private /*unused*/, double energy, double f, double m, Parity p,
                 std::string species, int n, double nu, double nui_exp, double nui_std,
                 double l_exp, double l_std, double s_exp, double s_std, double j_exp, double j_std,
                 double l_ryd_exp, double l_ryd_std, double j_ryd_exp, double j_ryd_std,
                 bool is_j_total_momentum, bool is_calculated_with_mqdt,
                 double underspecified_channel_contribution, Database &database,
                 size_t id_in_database)
    : Ket(energy, f, m, p), species(std::move(species)), quantum_number_n(n), quantum_number_nu(nu),
      quantum_number_nui_exp(nui_exp), quantum_number_nui_std(nui_std), quantum_number_l_exp(l_exp),
      quantum_number_l_std(l_std), quantum_number_s_exp(s_exp), quantum_number_s_std(s_std),
      quantum_number_j_exp(j_exp), quantum_number_j_std(j_std), quantum_number_l_ryd_exp(l_ryd_exp),
      quantum_number_l_ryd_std(l_ryd_std), quantum_number_j_ryd_exp(j_ryd_exp),
      quantum_number_j_ryd_std(j_ryd_std), is_j_total_momentum_(is_j_total_momentum),
      is_calculated_with_mqdt_(is_calculated_with_mqdt),
      underspecified_channel_contribution(underspecified_channel_contribution), database(database),
      id_in_database(id_in_database) {}

Database &KetAtom::get_database() const { return database; }

size_t KetAtom::get_id_in_database() const { return id_in_database; }

std::string KetAtom::get_label() const {
    size_t pos = species.find('_');
    std::string label = (pos != std::string::npos) ? species.substr(0, pos) : species;
    label[0] = static_cast<char>(std::toupper(label[0]));

    if (!is_calculated_with_mqdt_) {
        if (quantum_number_s_exp == 0) {
            label += "_singlet";
        } else if (quantum_number_s_exp == 1) {
            label += "_triplet";
        } else if (quantum_number_s_exp != 0.5) {
            throw std::runtime_error(
                "Invalid value for quantum number s in the single-channel description.");
        }
    }

    label += ":";

    if (is_calculated_with_mqdt_) {
        label += fmt::format("S={:.1f},nu={:.1f},L={:.1f},", quantum_number_s_exp,
                             quantum_number_nu, quantum_number_l_exp);
        label += this->is_j_total_momentum_ ? "J=" : "F=";
    } else {
        label += fmt::format("{:d},", quantum_number_n);
        if (quantum_number_l_exp == std::rint(quantum_number_l_exp) &&
            quantum_number_l_exp < quantum_number_l_labels.size()) {
            label += quantum_number_l_labels.at(static_cast<size_t>(quantum_number_l_exp));
        } else {
            label += fmt::format("{:.0f}", quantum_number_l_exp);
        }
        label += "_";
    }

    double total_momentum =
        this->is_j_total_momentum_ ? this->quantum_number_j_exp : this->quantum_number_f;
    if (total_momentum == std::rint(total_momentum)) {
        label += fmt::format("{:.0f}", total_momentum);
    } else if (2 * total_momentum == std::rint(2 * total_momentum)) {
        label += fmt::format("{:.0f}/2", 2 * total_momentum);
    } else {
        std::abort(); // can't happen because the total momentum is validated to be an integer
                      // or half-integer
    }

    if (this->quantum_number_m == std::rint(this->quantum_number_m)) {
        label += fmt::format(",{:.0f}", this->quantum_number_m);
    } else if (2 * this->quantum_number_m == std::rint(2 * this->quantum_number_m)) {
        label += fmt::format(",{:.0f}/2", 2 * this->quantum_number_m);
    } else {
        std::abort(); // can't happen because the quantum number m is validated to be an integer
                      // or half-integer
    }

    return label;
}

std::shared_ptr<KetAtom>
KetAtom::get_ket_for_different_quantum_number_m(double new_quantum_number_m) const {
    auto ket = *this;
    ket.quantum_number_m = new_quantum_number_m;
    return std::make_shared<KetAtom>(ket);
}

const std::string &KetAtom::get_species() const { return species; }

int KetAtom::get_quantum_number_n() const { return quantum_number_n; }

double KetAtom::get_quantum_number_nu() const { return quantum_number_nu; }

double KetAtom::get_quantum_number_nui() const { return quantum_number_nui_exp; }

double KetAtom::get_quantum_number_l() const { return quantum_number_l_exp; }

double KetAtom::get_quantum_number_s() const { return quantum_number_s_exp; }

double KetAtom::get_quantum_number_j() const { return quantum_number_j_exp; }

double KetAtom::get_quantum_number_l_ryd() const { return quantum_number_l_ryd_exp; }

double KetAtom::get_quantum_number_j_ryd() const { return quantum_number_j_ryd_exp; }

double KetAtom::get_quantum_number_nui_std() const { return quantum_number_nui_std; }

double KetAtom::get_quantum_number_l_std() const { return quantum_number_l_std; }

double KetAtom::get_quantum_number_s_std() const { return quantum_number_s_std; }

double KetAtom::get_quantum_number_j_std() const { return quantum_number_j_std; }

double KetAtom::get_quantum_number_l_ryd_std() const { return quantum_number_l_ryd_std; }

double KetAtom::get_quantum_number_j_ryd_std() const { return quantum_number_j_ryd_std; }

bool KetAtom::is_j_total_momentum() const { return is_j_total_momentum_; }

bool KetAtom::is_calculated_with_mqdt() const { return is_calculated_with_mqdt_; }

double KetAtom::get_underspecified_channel_contribution() const {
    return underspecified_channel_contribution;
}

bool KetAtom::operator==(const KetAtom &other) const {
    return Ket::operator==(other) && species == other.species &&
        quantum_number_n == other.quantum_number_n &&
        quantum_number_nu == other.quantum_number_nu &&
        quantum_number_nui_exp == other.quantum_number_nui_exp &&
        quantum_number_nui_std == other.quantum_number_nui_std &&
        quantum_number_l_exp == other.quantum_number_l_exp &&
        quantum_number_l_std == other.quantum_number_l_std &&
        quantum_number_s_exp == other.quantum_number_s_exp &&
        quantum_number_s_std == other.quantum_number_s_std &&
        quantum_number_j_exp == other.quantum_number_j_exp &&
        quantum_number_j_std == other.quantum_number_j_std;
}

bool KetAtom::operator!=(const KetAtom &other) const { return !(*this == other); }

size_t KetAtom::hash::operator()(const KetAtom &k) const {
    size_t seed = typename Ket::hash()(k);
    utils::hash_combine(seed, k.species);
    utils::hash_combine(seed, k.quantum_number_n);
    utils::hash_combine(seed, k.quantum_number_nu);
    utils::hash_combine(seed, k.quantum_number_nui_exp);
    utils::hash_combine(seed, k.quantum_number_nui_std);
    utils::hash_combine(seed, k.quantum_number_l_exp);
    utils::hash_combine(seed, k.quantum_number_l_std);
    utils::hash_combine(seed, k.quantum_number_s_exp);
    utils::hash_combine(seed, k.quantum_number_s_std);
    utils::hash_combine(seed, k.quantum_number_j_exp);
    utils::hash_combine(seed, k.quantum_number_j_std);
    return seed;
}

} // namespace pairinteraction
