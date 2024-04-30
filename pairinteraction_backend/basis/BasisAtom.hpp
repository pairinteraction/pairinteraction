#pragma once

#include "basis/Basis.hpp"
#include "utils/traits.hpp"

#include <complex>
#include <string>
#include <vector>

class Database;

template <typename T>
class BasisAtomCreator;

template <typename Real>
class KetAtom;

// Specialize BasisTraits for BasisAtom
template <typename Scalar>
class BasisAtom;

template <typename Scalar>
struct traits::BasisTraits<BasisAtom<Scalar>> {
    using scalar_t = Scalar;
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetAtom<real_t>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;
};

/**
 * @class BasisAtom
 *
 * @brief Class for creating a basis of atomic kets.
 *
 * @tparam Scalar Complex number type.
 */
template <typename Scalar>
class BasisAtom : public Basis<BasisAtom<Scalar>> {
public:
    using Type = BasisAtom<Scalar>;
    using ketvec_t = typename traits::BasisTraits<Type>::ketvec_t;

    Database &get_database() const;

private:
    friend class Database;
    BasisAtom(ketvec_t &&kets, std::string table, Database &database, std::string species);
    std::string table;
    Database &database;
    std::string species;
};

extern template class BasisAtom<float>;
extern template class BasisAtom<double>;
extern template class BasisAtom<std::complex<float>>;
extern template class BasisAtom<std::complex<double>>;
