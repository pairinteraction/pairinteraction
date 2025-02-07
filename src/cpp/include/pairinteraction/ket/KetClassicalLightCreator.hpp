#pragma once

#include <memory>
#include <optional>
#include <type_traits>

namespace pairinteraction {
class KetClassicalLight;

/**
 * @class KetClassicalLight
 *
 * @brief Builder class for creating KetClassicalLight.
 */
class KetClassicalLightCreator {
public:
    KetClassicalLightCreator() = default;
    KetClassicalLightCreator(double energy, int q);
    KetClassicalLightCreator &set_photon_energy(double value);
    KetClassicalLightCreator &set_quantum_number_q(int value);
    std::shared_ptr<const KetClassicalLight> create() const;

private:
    std::optional<double> photon_energy;
    std::optional<int> quantum_number_q;
};
} // namespace pairinteraction
