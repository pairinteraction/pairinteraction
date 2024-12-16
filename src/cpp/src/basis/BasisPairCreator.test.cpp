#include "pairinteraction/basis/BasisPairCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerEigen.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/streamed.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("create a BasisPair") {
    // Create single-atom system
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 2)
                     .create(database);
    SystemAtom<double> system(basis);
    system.set_electric_field({0, 0, 1 * 1.9446903811524456e-10});

    DiagonalizerEigen<double> diagonalizer;
    system.diagonalize(diagonalizer);

    // Get energy window for a two-atom basis
    auto ket = KetAtomCreator<double>()
                   .set_species("Rb")
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(0)
                   .set_quantum_number_m(0.5)
                   .create(database);
    double min_energy = 2 * ket->get_energy() - 3 / 6579683.920501762;
    double max_energy = 2 * ket->get_energy() + 3 / 6579683.920501762;

    // Create two-atom bases
    auto basis_pair_a = pairinteraction::BasisPairCreator<double>()
                            .add(system)
                            .add(system)
                            .restrict_energy(min_energy, max_energy)
                            .restrict_quantum_number_m(1, 1)
                            .create();
    auto basis_pair_b = pairinteraction::BasisPairCreator<double>()
                            .add(system)
                            .add(system)
                            .restrict_energy(min_energy, max_energy)
                            .restrict_quantum_number_m(1, 1)
                            .create();

    DOCTEST_SUBCASE("check equality of kets") {
        // Obtain kets from the two-atom bases and check for equality
        auto ket1a = basis_pair_a->get_kets()[0];
        auto ket1b = basis_pair_b->get_kets()[0];
        auto ket2a = basis_pair_a->get_kets()[1];
        auto ket2b = basis_pair_b->get_kets()[1];
        DOCTEST_CHECK(*ket1a == *ket1a);
        DOCTEST_CHECK(*ket2a == *ket2a);
        DOCTEST_CHECK(*ket1a != *ket2b);
        DOCTEST_CHECK(*ket2a != *ket1b);

        // Currently, kets from different BasisPair are never equal
        DOCTEST_CHECK(*ket1a != *ket1b);
        DOCTEST_CHECK(*ket2a != *ket2b);
    }

    DOCTEST_SUBCASE("check overlap") {
        auto overlaps = basis_pair_a->get_overlaps(ket, ket);

        // The total overlap is less than 1 because of the restricted energy window
        DOCTEST_CHECK(overlaps.sum() == doctest::Approx(0.9107819201));
    }

    DOCTEST_SUBCASE("get the atomic states constituting a ket of the basis_pair") {
        auto atomic_states = basis_pair_a->get_kets()[0]->get_atomic_states();
        DOCTEST_CHECK(atomic_states.size() == 2);
        DOCTEST_CHECK(atomic_states[0]->get_number_of_states() == 1);
        DOCTEST_CHECK(atomic_states[0]->get_number_of_kets() == basis->get_number_of_kets());
    }
}
} // namespace pairinteraction
