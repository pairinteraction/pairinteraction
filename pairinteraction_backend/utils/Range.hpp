#pragma once

#include "utils/traits.hpp"

template <typename Sortable>
class Range {
    static_assert(traits::OpTraits<Sortable>::has_less_v);

public:
    Range() = default;
    Range(Sortable min, Sortable max) : min_(min), max_(max), is_finite_(true) {
        if (max < min) {
            throw std::invalid_argument("It must be max >= min.");
        }
    }
    Sortable min() const {
        if (!is_finite_) {
            throw std::runtime_error("The range is infinite.");
        }
        return min_;
    }
    Sortable max() const {
        if (!is_finite_) {
            throw std::runtime_error("The range is infinite.");
        }
        return max_;
    }
    bool is_finite() const { return is_finite_; }

private:
    Sortable min_{};
    Sortable max_{};
    bool is_finite_{false};
};
