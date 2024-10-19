#include "pairinteraction/basis/BasisAtom.hpp"

#include "pairinteraction/ket/KetAtom.hpp"

namespace pairinteraction {
template <typename Scalar>
BasisAtom<Scalar>::BasisAtom(Private /*unused*/, ketvec_t &&kets, std::string &&id_of_kets,
                             Database &database)
    : Basis<BasisAtom<Scalar>>(std::move(kets), std::move(id_of_kets)), database(database) {}

template <typename Scalar>
Database &BasisAtom<Scalar>::get_database() const {
    return database;
}

template <typename Scalar>
const std::string &BasisAtom<Scalar>::get_species() const {
    return this->kets[0]->get_species();
}

// Explicit instantiations
template class BasisAtom<float>;
template class BasisAtom<double>;
template class BasisAtom<std::complex<float>>;
template class BasisAtom<std::complex<double>>;
} // namespace pairinteraction
