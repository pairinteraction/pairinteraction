#pragma once

#include "pairinteraction/ket/Ket.hpp"

#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace pairinteraction {
class Database;

template <typename Real>
class KetAtom;

template <typename Real>
class KetCombined : public Ket<Real> {
    static_assert(std::is_floating_point_v<Real>);

public:
    // kets.emplace_back(id, energy, quantum_number_f, quantum_number_m, parity, {ket1, ket2});

    KetCombined(size_t id, Real energy, Real quantum_number_f, Real quantum_number_m, Parity parity,
                std::vector<std::shared_ptr<const KetAtom<Real>>> &&kets_with_largest_overlap);
    std::string get_label() const override;
    size_t get_id() const override;
    size_t get_id_for_different_quantum_number_m(Real new_quantum_number_m) const override;

private:
    size_t id;
    std::vector<std::shared_ptr<const KetAtom<Real>>> kets_with_largest_overlap;
};

extern template class KetCombined<float>;
extern template class KetCombined<double>;
} // namespace pairinteraction
