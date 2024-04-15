#include "ket/Ket.hpp"

#include <limits>

template <typename Real>
Ket<Real>::Ket(Real energy, float f, float m, int p, std::string label)
    : energy(energy), quantum_number_f(f), quantum_number_m(m), parity(p), label(label) {}

template <typename Real>
Real Ket<Real>::get_energy() const {
    return energy;
}

template <typename Real>
float Ket<Real>::get_quantum_number_f() const {
    return quantum_number_f;
}

template <typename Real>
float Ket<Real>::get_quantum_number_m() const {
    return quantum_number_m;
}

template <typename Real>
int Ket<Real>::get_parity() const {
    return parity;
}

template <typename Real>
std::string Ket<Real>::get_label() const {
    return label;
}

// Explicit instantiations
template class Ket<float>;
template class Ket<double>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>
#include <sstream>

DOCTEST_TEST_CASE("constructing a class derived from ket") {
    class KetDerived : public Ket<float> {
    private:
        friend class KetDerivedCreator;
        KetDerived(float energy, float f, float m, int p, std::string label)
            : Ket<float>(energy, f, m, p, label) {}
    };

    class KetDerivedCreator {
    public:
        KetDerivedCreator(float energy, float f, float m, int p, std::string label)
            : energy(energy), f(f), m(m), p(p), label(label) {}
        KetDerived create() const { return KetDerived(energy, f, m, p, label); }

    private:
        float energy;
        float f;
        float m;
        int p;
        std::string label;
    };

    auto ket = KetDerivedCreator(1.0f, 2.0f, 3.0f, 4, "my_label").create();

    // Check that the label can be printed
    std::stringstream ss;
    ss << ket;
    DOCTEST_CHECK(ss.str() == "my_label");

    // Output the label to the doctest log
    //SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", ket); // FIXME: Cannot format an argument.
}
