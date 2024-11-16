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
 *
 * @tparam Real Real number type.
 */
template <typename Real>
class KetAtom : public Ket<Real> {
    static_assert(std::is_floating_point_v<Real>);

    friend class Database;
    struct Private {};

public:
    KetAtom(Private /*unused*/, Real energy, Real f, Real m, Parity p, std::string species, int n,
            Real nu_exp, Real nu_std, Real l_exp, Real l_std, Real s_exp, Real s_std, Real j_exp,
            Real j_std, Database &database, size_t id_in_database);

    Database &get_database() const;
    size_t get_id_in_database() const;
    std::string get_label() const override;
    std::shared_ptr<KetAtom<Real>>
    get_ket_for_different_quantum_number_m(Real new_quantum_number_m) const;
    const std::string &get_species() const;
    int get_quantum_number_n() const;
    Real get_quantum_number_nu() const;
    Real get_quantum_number_l() const;
    Real get_quantum_number_s() const;
    Real get_quantum_number_j() const;
    Real get_quantum_number_nu_std() const;
    Real get_quantum_number_l_std() const;
    Real get_quantum_number_s_std() const;
    Real get_quantum_number_j_std() const;

    bool operator==(const KetAtom<Real> &other) const;

    struct hash {
        std::size_t operator()(const KetAtom<Real> &k) const;
    };

private:
    std::string species;
    int quantum_number_n;
    Real quantum_number_nu_exp;
    Real quantum_number_nu_std;
    Real quantum_number_l_exp;
    Real quantum_number_l_std;
    Real quantum_number_s_exp;
    Real quantum_number_s_std;
    Real quantum_number_j_exp;
    Real quantum_number_j_std;
    Database &database;
    size_t id_in_database;
};

extern template class KetAtom<float>;
extern template class KetAtom<double>;
} // namespace pairinteraction
