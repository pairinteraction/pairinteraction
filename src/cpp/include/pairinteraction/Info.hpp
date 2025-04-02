// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

namespace pairinteraction {

struct Info {
    // Eigen diagonalizer is always available
    bool has_eigen = true;

    // LAPACKE diagonalizers are available if compiled with either MKL or LAPACKE
#if defined(WITH_MKL) || defined(WITH_LAPACKE)
    bool has_lapacke_evd = true;
    bool has_lapacke_evr = true;
#else
    bool has_lapacke_evd = false;
    bool has_lapacke_evr = false;
#endif

    // FEAST diagonalizer is only available if compiled with MKL
#if defined(WITH_MKL)
    bool has_feast = true;
#else
    bool has_feast = false;
#endif
};

} // namespace pairinteraction
