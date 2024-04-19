#pragma once

#include <string>

#include "ket/Ket.hpp"

class Database;

template <typename Real>
class KetAtomCreator;

/**
 * @class KetAtom
 *
 * @brief Class for representing atomic kets.
 *
 * @tparam Real Real number type.
 */
template <typename Real>
class KetAtom : public Ket<Real> {
public:
    std::string get_species() const;
    int get_quantum_number_n() const;
    Real get_quantum_number_nu() const;
    Real get_quantum_number_l() const;
    Real get_quantum_number_s() const;
    Real get_quantum_number_j() const;

private:
    friend class KetAtomCreator<Real>;
    KetAtom(Real energy, float f, float m, int p, std::string label, size_t id, std::string species,
            int n, Real nu_exp, Real nu_std, Real l_exp, Real l_std, Real s_exp, Real s_std,
            Real j_exp, Real j_std, Database &database);
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
};

extern template class KetAtom<float>;
extern template class KetAtom<double>;
