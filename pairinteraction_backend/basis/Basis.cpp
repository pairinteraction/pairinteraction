#include "Basis.hpp"

Basis::Basis()
{
}

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>

DOCTEST_TEST_CASE("constructing a class derived from basis") {
    class DerivedBasis : public Basis
    {
    public:
        DerivedBasis() : Basis() {}
    };

    DerivedBasis derivedBasis;
    DOCTEST_CHECK(false);
}

DOCTEST_TEST_CASE("constructing a class derived from basis2") {
    class DerivedBasis : public Basis
    {
    public:
        DerivedBasis() : Basis() {}
    };

    DerivedBasis derivedBasis;
    DOCTEST_CHECK(false);
}