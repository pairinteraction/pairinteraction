// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <cpptrace/cpptrace.hpp>

#ifdef eigen_assert
#error "'utils/eigen_assertion.hpp' must be included before any Eigen header"
#else
#define eigen_assert(X)                                                                            \
    if (!(X)) {                                                                                    \
        throw cpptrace::runtime_error("<EIGEN> Assertion " EIGEN_MAKESTRING(                       \
            X) " failed at " __FILE__ ":" EIGEN_MAKESTRING(__LINE__));                             \
    }
#endif
