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
