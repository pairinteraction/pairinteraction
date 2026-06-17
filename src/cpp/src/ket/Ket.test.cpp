// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/Ket.hpp"

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
        KetDerived(Private /*unused*/, double m) : Ket(0), quantum_number_m(m) {}
        std::string get_label() const override { return "my_label"; }
        std::shared_ptr<KetDerived>
        get_ket_for_different_quantum_number_m(double new_quantum_number_m) const {
            auto ket = *this;
            ket.quantum_number_m = new_quantum_number_m;
            return std::make_shared<KetDerived>(ket);
        }

    private:
        double quantum_number_m;
    };

    class KetDerivedCreator {
    public:
        KetDerivedCreator(double m) : m(m) {}
        std::shared_ptr<const KetDerived> create() const {
            return std::make_shared<const KetDerived>(KetDerived::Private(), m);
        }

    private:
        double m;
    };

    auto ket = KetDerivedCreator(3.0).create();

    // Check that the label can be printed
    std::stringstream ss;
    ss << *ket;
    DOCTEST_CHECK(ss.str() == "my_label");

    // Output the label to the doctest log
    DOCTEST_MESSAGE("Ket: ", *ket);
}
} // namespace pairinteraction
