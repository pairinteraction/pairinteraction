#include "ket/Ket.hpp"

#include <limits>

template <typename Real>
Ket<Real>::Ket(Real energy, Real f, Real m, int p)
    : energy(energy), quantum_number_f(f), quantum_number_m(m), parity(p) {}

template <typename Real>
Real Ket<Real>::get_energy() const {
    return energy;
}

template <typename Real>
Real Ket<Real>::get_quantum_number_f() const {
    return quantum_number_f;
}

template <typename Real>
Real Ket<Real>::get_quantum_number_m() const {
    return quantum_number_m;
}

template <typename Real>
int Ket<Real>::get_parity() const {
    return parity;
}

// Explicit instantiations
template class Ket<float>;
template class Ket<double>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include "utils/hash.hpp"
#include "utils/streamed.hpp"
#include <doctest/doctest.h>
#include <spdlog/spdlog.h>
#include <sstream>

DOCTEST_TEST_CASE("constructing a class derived from ket") {
    class KetDerivedCreator;

    class KetDerived : public Ket<float> {
        friend class KetDerivedCreator;
        struct Private {};

    public:
        KetDerived(Private /*unused*/, float f, float m, int p) : Ket<float>(0, f, m, p) {}
        std::string get_label() const override { return "my_label"; }
        size_t get_id() const override {
            size_t seed = 0;
            hash::hash_combine(seed, this->quantum_number_f);
            hash::hash_combine(seed, this->quantum_number_m);
            hash::hash_combine(seed, this->parity);
            return seed;
        }
        size_t get_id_for_different_quantum_number_m(float new_quantum_number_m) const override {
            size_t seed = 0;
            hash::hash_combine(seed, this->quantum_number_f);
            hash::hash_combine(seed, new_quantum_number_m);
            hash::hash_combine(seed, this->parity);
            return seed;
        }
    };

    class KetDerivedCreator {
    public:
        KetDerivedCreator(float f, float m, int p) : f(f), m(m), p(p) {}
        std::shared_ptr<const KetDerived> create() const {
            return std::make_shared<const KetDerived>(KetDerived::Private(), f, m, p);
        }

    private:
        float f;
        float m;
        int p;
    };

    auto ket = KetDerivedCreator(2.0F, 3.0F, 4).create();

    // Check that the label can be printed
    std::stringstream ss;
    ss << *ket;
    DOCTEST_CHECK(ss.str() == "my_label");

    // Output the label to the doctest log
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(*ket));
}
