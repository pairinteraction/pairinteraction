// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

namespace pairinteraction {

struct Info {
    Info() = delete;

    // Whether it compiled with MKL
#if defined(WITH_MKL)
    static constexpr bool with_mkl = true;
#else
    static constexpr bool with_mkl = false;
#endif

    // Eigen diagonalizer is always available
    static constexpr bool has_eigen = true;

    // LAPACKE diagonalizers are available if compiled with either MKL or LAPACKE
#if defined(WITH_MKL) || defined(WITH_LAPACKE)
    static constexpr bool has_lapacke_evd = true;
    static constexpr bool has_lapacke_evr = true;
#else
    static constexpr bool has_lapacke_evd = false;
    static constexpr bool has_lapacke_evr = false;
#endif

    // FEAST diagonalizer is only available if compiled with MKL
#if defined(WITH_MKL)
    static constexpr bool has_feast = true;
#else
    static constexpr bool has_feast = false;
#endif
};

} // namespace pairinteraction
