#pragma once

#include <exception>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>

namespace pintr {
template <typename Real>
class KetClassicalLight;

/**
 * @class KetClassicalLight
 *
 * @brief Builder class for creating KetClassicalLight.
 *
 * @tparam Real Real number type.
 */
template <typename Real>
class KetClassicalLightCreator {
    static_assert(std::is_floating_point_v<Real>);

public:
    KetClassicalLightCreator() = default;
    KetClassicalLightCreator(Real energy, int q);
    KetClassicalLightCreator<Real> &set_photon_energy(Real value);
    KetClassicalLightCreator<Real> &set_quantum_number_q(int value);
    std::shared_ptr<const KetClassicalLight<Real>> create() const;

private:
    std::optional<Real> photon_energy;
    std::optional<int> quantum_number_q;
};

extern template class KetClassicalLightCreator<float>;
extern template class KetClassicalLightCreator<double>;
} // namespace pintr
