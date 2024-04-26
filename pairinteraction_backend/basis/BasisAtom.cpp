#include "basis/BasisAtom.hpp"

template <typename Scalar>
BasisAtom<Scalar>::BasisAtom(ketvec_t &&kets, std::string table, Database &database,
                             std::string species)
    : Basis<BasisAtom<Scalar>>(std::move(kets)), table(table), database(database),
      species(species) {}

// Explicit instantiations
template class BasisAtom<float>;
template class BasisAtom<double>;
template class BasisAtom<std::complex<float>>;
template class BasisAtom<std::complex<double>>;
