#pragma once

#include <exception>
#include <limits>
#include <optional>
#include <string>

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
public:
    KetClassicalLightCreator() = default;
    KetClassicalLightCreator(Real energy, int q);
    KetClassicalLightCreator<Real> &set_energy(Real value);
    KetClassicalLightCreator<Real> &set_quantum_number_q(int value);
    KetClassicalLight<Real> create() const;

private:
    std::optional<Real> energy;
    std::optional<int> quantum_number_q;
};

extern template class KetClassicalLightCreator<float>;
extern template class KetClassicalLightCreator<double>;
