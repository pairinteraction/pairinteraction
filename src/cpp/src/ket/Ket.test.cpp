#include "pairinteraction/ket/Ket.hpp"

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <doctest/doctest.h>
#include <memory>
#include <sstream>

namespace pairinteraction {
DOCTEST_TEST_CASE("constructing a class derived from ket") {
    class KetDerivedCreator;

    class KetDerived : public Ket<float> {
        friend class KetDerivedCreator;
        struct Private {};

    public:
        KetDerived(Private /*unused*/, float f, float m, Parity p) : Ket<float>(0, f, m, p) {}
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
        KetDerivedCreator(float f, float m, Parity p) : f(f), m(m), p(p) {}
        std::shared_ptr<const KetDerived> create() const {
            return std::make_shared<const KetDerived>(KetDerived::Private(), f, m, p);
        }

    private:
        float f;
        float m;
        Parity p;
    };

    auto ket = KetDerivedCreator(2.0F, 3.0F, Parity::EVEN).create();

    // Check that the label can be printed
    std::stringstream ss;
    ss << *ket;
    DOCTEST_CHECK(ss.str() == "my_label");

    // Output the label to the doctest log
    DOCTEST_MESSAGE("Ket: ", *ket);
}
} // namespace pairinteraction
