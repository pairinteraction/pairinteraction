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
                real_t energy, real_t quantum_number_f, real_t quantum_number_m, Parity parity);
    std::string get_label() const override;
    size_t get_id() const override;
    size_t get_id_for_different_quantum_number_m(real_t new_quantum_number_m) const override;

private:
    std::vector<size_t> atomic_indices;
    std::vector<std::string> atomic_labels;
    std::vector<std::shared_ptr<const BasisAtom<Scalar>>> atomic_bases;
};

extern template class KetCombined<float>;
extern template class KetCombined<double>;
extern template class KetCombined<std::complex<float>>;
extern template class KetCombined<std::complex<double>>;
} // namespace pairinteraction
