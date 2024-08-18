#include "pairinteraction/basis/BasisAtom.hpp"

namespace pairinteraction {
template <typename Scalar>
BasisAtom<Scalar>::BasisAtom(Private /*unused*/, ketvec_t &&kets, std::string table,
                             Database &database, std::string species)
    : Basis<BasisAtom<Scalar>>(std::move(kets)), table(std::move(table)), database(database),
      species(std::move(species)) {}

template <typename Scalar>
Database &BasisAtom<Scalar>::get_database() const {
    return database;
}

// Explicit instantiations
template class BasisAtom<float>;
template class BasisAtom<double>;
template class BasisAtom<std::complex<float>>;
template class BasisAtom<std::complex<double>>;
} // namespace pairinteraction
