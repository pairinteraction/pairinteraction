// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/Ket.hpp"

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <doctest/doctest.h>
#include <memory>
#include <sstream>

namespace pairinteraction {
DOCTEST_TEST_CASE("constructing a class derived from ket") {
    class KetDerivedCreator;

    class KetDerived : public Ket {
        friend class KetDerivedCreator;
        struct Private {};

    public:
        KetDerived(Private /*unused*/, double f, double m, Parity p) : Ket(0, f, m, p) {}
        std::string get_label() const override { return "my_label"; }
        std::shared_ptr<KetDerived>
        get_ket_for_different_quantum_number_m(double new_quantum_number_m) const {
            auto ket = *this;
            ket.quantum_number_m = new_quantum_number_m;
            return std::make_shared<KetDerived>(ket);
        }
    };

    class KetDerivedCreator {
    public:
        KetDerivedCreator(double f, double m, Parity p) : f(f), m(m), p(p) {}
        std::shared_ptr<const KetDerived> create() const {
            return std::make_shared<const KetDerived>(KetDerived::Private(), f, m, p);
        }

    private:
        double f;
        double m;
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
