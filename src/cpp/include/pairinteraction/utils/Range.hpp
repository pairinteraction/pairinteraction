// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/traits.hpp"

namespace pairinteraction {
template <typename Sortable>
class Range {
    static_assert(traits::OpTraits<Sortable>::has_less_v);

public:
    Range() = default;
    Range(Sortable min, Sortable max) : _min(min), _max(max), _is_finite(true) {
        if (max < min) {
            throw std::invalid_argument("It must be max >= min.");
        }
    }
    Sortable min() const {
        if (!_is_finite) {
            throw std::runtime_error("The range is infinite.");
        }
        return _min;
    }
    Sortable max() const {
        if (!_is_finite) {
            throw std::runtime_error("The range is infinite.");
        }
        return _max;
    }
    bool is_finite() const { return _is_finite; }

private:
    Sortable _min{};
    Sortable _max{};
    bool _is_finite{false};
};
} // namespace pairinteraction
