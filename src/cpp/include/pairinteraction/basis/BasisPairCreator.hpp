// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <memory>
#include <optional>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class BasisPair;

template <typename Scalar>
class BasisAtom;

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
    BasisPairCreator<Scalar> &add(std::shared_ptr<const BasisAtom<Scalar>> basis_atom);
    BasisPairCreator<Scalar> &restrict_energy(real_t min, real_t max);
    BasisPairCreator<Scalar> &restrict_quantum_number_m(real_t min, real_t max);
    BasisPairCreator<Scalar> &restrict_parity_under_inversion(int value);
    BasisPairCreator<Scalar> &restrict_parity_under_permutation(int value);
    std::shared_ptr<const BasisPair<Scalar>> create() const;

private:
    std::vector<std::shared_ptr<const BasisAtom<Scalar>>> bases_atom;
    Range<real_t> range_energy;
    Range<real_t> range_quantum_number_m;
    std::optional<int> parity_under_inversion;
    std::optional<int> parity_under_permutation;
};

extern template class BasisPairCreator<double>;
extern template class BasisPairCreator<std::complex<double>>;
} // namespace pairinteraction
