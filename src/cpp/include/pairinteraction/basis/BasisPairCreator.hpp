// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <functional>
#include <memory>
#include <vector>

namespace pairinteraction {
enum class Parity : int;

template <typename Scalar>
class BasisPair;

template <typename Scalar>
class SystemAtom;

template <typename Scalar>
class KetPair;

template <typename Scalar>
class BasisPairCreator {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

public:
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using basis_t = BasisPair<Scalar>;
    using ket_t = KetPair<Scalar>;
    using ketvec_t = std::vector<std::shared_ptr<const ket_t>>;

    BasisPairCreator() = default;
    BasisPairCreator<Scalar> &add(const SystemAtom<Scalar> &system_atom);
    BasisPairCreator<Scalar> &restrict_energy(real_t min, real_t max);
    BasisPairCreator<Scalar> &restrict_quantum_number_m(real_t min, real_t max);
    BasisPairCreator<Scalar> &restrict_parity_under_inversion(Parity value);
    BasisPairCreator<Scalar> &restrict_parity_under_permutation(Parity value);
    std::shared_ptr<const BasisPair<Scalar>> create() const;

private:
    std::vector<std::reference_wrapper<const SystemAtom<Scalar>>> systems_atom;
    Range<real_t> range_energy;
    Range<real_t> range_quantum_number_m;
    Parity parity_under_inversion{Parity::UNKNOWN};
    Parity parity_under_permutation{Parity::UNKNOWN};
};

extern template class BasisPairCreator<double>;
extern template class BasisPairCreator<std::complex<double>>;
} // namespace pairinteraction
