#pragma once

namespace pairinteraction {

struct Info {
    // Eigen diagonalizer is always available
    static inline constexpr bool has_eigen = true;

    // LAPACKE diagonalizers are available if compiled with either MKL or LAPACKE
#if defined(WITH_MKL) || defined(WITH_LAPACKE)
    static inline constexpr bool has_lapacke_evd = true;
    static inline constexpr bool has_lapacke_evr = true;
#else
    static inline constexpr bool has_lapacke_evd = false;
    static inline constexpr bool has_lapacke_evr = false;
#endif

    // FEAST diagonalizer is only available if compiled with MKL
#if defined(WITH_MKL)
    static inline constexpr bool has_feast = true;
#else
    static inline constexpr bool has_feast = false;
#endif
};

} // namespace pairinteraction
