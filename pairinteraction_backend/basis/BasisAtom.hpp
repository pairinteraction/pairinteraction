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

// Specialize CrtpTraits for BasisAtom
template <typename Scalar>
class BasisAtom;

template <typename Scalar>
struct traits::CrtpTraits<BasisAtom<Scalar>> {
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
    static_assert(traits::is_complex_or_floating_point_v<Scalar>);

    friend class Database;
    struct Private {};

public:
    using Type = BasisAtom<Scalar>;
    using ketvec_t = typename traits::CrtpTraits<Type>::ketvec_t;

    BasisAtom(Private, ketvec_t &&kets, std::string table, Database &database, std::string species);
    Database &get_database() const;

private:
    std::string table;
    Database &database;
    std::string species;
};

extern template class BasisAtom<float>;
extern template class BasisAtom<double>;
extern template class BasisAtom<std::complex<float>>;
extern template class BasisAtom<std::complex<double>>;
