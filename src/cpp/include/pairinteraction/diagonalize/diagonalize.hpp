// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <functional>
#include <initializer_list>
#include <optional>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class DiagonalizerInterface;

// Note that although a vector is also constructible from a list, the overload resolution
// will prefer the initializer list overload because of less conversions required.

template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 std::optional<typename Derived::real_t> min_eigenenergy = {},
                 std::optional<typename Derived::real_t> max_eigenenergy = {}, double rtol = 1e-6);

template <typename Derived>
void diagonalize(std::vector<Derived> &systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 std::optional<typename Derived::real_t> min_eigenenergy = {},
                 std::optional<typename Derived::real_t> max_eigenenergy = {}, double rtol = 1e-6);

template <typename Derived>
void diagonalize(std::vector<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 std::optional<typename Derived::real_t> min_eigenenergy = {},
                 std::optional<typename Derived::real_t> max_eigenenergy = {}, double rtol = 1e-6);

} // namespace pairinteraction
