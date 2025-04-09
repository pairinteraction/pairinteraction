// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/ket/Ket.hpp"

#include <string>
#include <type_traits>

namespace pairinteraction {
class Database;

enum class Parity : int;

/**
 * @class KetAtom
 *
 * @brief Class for representing atomic kets.
 */
class KetAtom : public Ket {
    friend class Database;
    struct Private {};

public:
    KetAtom(Private /*unused*/, double energy, double f, double m, Parity p, std::string species,
            int n, double nu, double nui_exp, double nui_std, double l_exp, double l_std,
            double s_exp, double s_std, double j_exp, double j_std, double l_ryd_exp,
            double l_ryd_std, double j_ryd_exp, double j_ryd_std, bool is_j_total_momentum,
            bool is_calculated_with_mqdt, double underspecified_channel_contribution,
            Database &database, size_t id_in_database);

    Database &get_database() const;
    size_t get_id_in_database() const;
    std::string get_label() const override;
    std::shared_ptr<KetAtom>
    get_ket_for_different_quantum_number_m(double new_quantum_number_m) const;
    const std::string &get_species() const;
    int get_quantum_number_n() const;
    double get_quantum_number_nu() const;
    double get_quantum_number_nui() const;
    double get_quantum_number_l() const;
    double get_quantum_number_s() const;
    double get_quantum_number_j() const;
    double get_quantum_number_l_ryd() const;
    double get_quantum_number_j_ryd() const;
    double get_quantum_number_nui_std() const;
    double get_quantum_number_l_std() const;
    double get_quantum_number_s_std() const;
    double get_quantum_number_j_std() const;
    double get_quantum_number_l_ryd_std() const;
    double get_quantum_number_j_ryd_std() const;
    bool is_j_total_momentum() const;
    bool is_calculated_with_mqdt() const;
    double get_underspecified_channel_contribution() const;

    bool operator==(const KetAtom &other) const;
    bool operator!=(const KetAtom &other) const;

    struct hash {
        std::size_t operator()(const KetAtom &k) const;
    };

private:
    std::string species;
    int quantum_number_n;
    double quantum_number_nu;
    double quantum_number_nui_exp;
    double quantum_number_nui_std;
    double quantum_number_l_exp;
    double quantum_number_l_std;
    double quantum_number_s_exp;
    double quantum_number_s_std;
    double quantum_number_j_exp;
    double quantum_number_j_std;
    double quantum_number_l_ryd_exp;
    double quantum_number_l_ryd_std;
    double quantum_number_j_ryd_exp;
    double quantum_number_j_ryd_std;
    bool is_j_total_momentum_;
    bool is_calculated_with_mqdt_;
    double underspecified_channel_contribution;
    Database &database;
    size_t id_in_database;
};

} // namespace pairinteraction
