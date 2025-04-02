// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <fmt/ostream.h>
#include <utility>

#if FMT_VERSION < 90000

namespace fmt {
template <typename T>
inline auto streamed(T &&v) {
    return std::forward<T>(v);
}
} // namespace fmt
#endif
