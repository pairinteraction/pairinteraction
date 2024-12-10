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
class BasisCombinedCreator;

template <typename Scalar>
class BasisAtom;

template <typename Scalar>
class KetCombined : public Ket<typename traits::NumTraits<Scalar>::real_t> {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using real_t = typename traits::NumTraits<Scalar>::real_t;

    friend class BasisCombinedCreator<Scalar>;
    struct Private {};

public:
    KetCombined(Private /*unused*/, std::initializer_list<size_t> atomic_indices,
                std::initializer_list<std::string> atomic_labels,
                std::initializer_list<std::shared_ptr<const BasisAtom<Scalar>>> atomic_bases,
                real_t energy);

    std::string get_label() const override;
    std::shared_ptr<KetCombined<Scalar>>
    get_ket_for_different_quantum_number_m(real_t new_quantum_number_m) const;
    std::vector<std::shared_ptr<const BasisAtom<Scalar>>> get_atomic_states() const;

    bool operator==(const KetCombined<Scalar> &other) const;
    bool operator!=(const KetCombined<Scalar> &other) const;

    struct hash {
        std::size_t operator()(const KetCombined<Scalar> &k) const;
    };

private:
    std::vector<size_t> atomic_indices;
    std::vector<std::string> atomic_labels;
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

extern template class KetCombined<float>;
extern template class KetCombined<double>;
extern template class KetCombined<std::complex<float>>;
extern template class KetCombined<std::complex<double>>;
} // namespace pairinteraction
