#pragma once

#include <stdexcept>

#ifdef eigen_assert
#error "'utils/eigen_assertion.hpp' must be included before any Eigen header"
#else
#define eigen_assert(X)                                                                            \
    if (!(X)) {                                                                                    \
        throw std::runtime_error("<EIGEN> Assertion " EIGEN_MAKESTRING(                            \
            X) " failed at " __FILE__ ":" EIGEN_MAKESTRING(__LINE__));                             \
    }
#endif
