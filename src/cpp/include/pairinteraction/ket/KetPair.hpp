// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/ket/Ket.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <complex>
#include <initializer_list>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class BasisPairCreator;

template <typename Scalar>
class BasisAtom;

template <typename Scalar>
class KetPair : public Ket {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using real_t = typename traits::NumTraits<Scalar>::real_t;

    friend class BasisPairCreator<Scalar>;
    struct Private {};

public:
    KetPair(Private /*unused*/, std::initializer_list<size_t> atomic_indices,
            std::initializer_list<std::shared_ptr<const BasisAtom<Scalar>>> atomic_bases,
            real_t energy);

    std::string get_label() const override;
    std::shared_ptr<KetPair<Scalar>>
    get_ket_for_different_quantum_number_m(real_t new_quantum_number_m) const;
    std::vector<std::shared_ptr<const BasisAtom<Scalar>>> get_atomic_states() const;

    bool operator==(const KetPair<Scalar> &other) const;
    bool operator!=(const KetPair<Scalar> &other) const;

    struct hash {
        std::size_t operator()(const KetPair<Scalar> &k) const;
    };

private:
    std::vector<size_t> atomic_indices;
    std::vector<std::shared_ptr<const BasisAtom<Scalar>>> atomic_bases;
    static real_t
    calculate_quantum_number_f(const std::vector<size_t> &indices,
                               const std::vector<std::shared_ptr<const BasisAtom<Scalar>>> &bases);
    static real_t
    calculate_quantum_number_m(const std::vector<size_t> &indices,
                               const std::vector<std::shared_ptr<const BasisAtom<Scalar>>> &bases);
    static Parity
    calculate_parity(const std::vector<size_t> &indices,
                     const std::vector<std::shared_ptr<const BasisAtom<Scalar>>> &bases);
};

extern template class KetPair<double>;
extern template class KetPair<std::complex<double>>;
} // namespace pairinteraction
