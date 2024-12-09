#include "pairinteraction/basis/BasisCombinedCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/basis/BasisCombined.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/ket/KetCombined.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/streamed.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("create a combined basis and check equality of kets") {
    // Create single-atom system
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 2)
                     .create(database);
    SystemAtom<double> system(basis);

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
    auto combined_basis_a = pairinteraction::BasisCombinedCreator<double>()
                                .add(system)
                                .add(system)
                                .restrict_energy(min_energy, max_energy)
                                .restrict_quantum_number_m(1, 1)
                                .create();
    auto combined_basis_b = pairinteraction::BasisCombinedCreator<double>()
                                .add(system)
                                .add(system)
                                .restrict_energy(min_energy, max_energy)
                                .restrict_quantum_number_m(1, 1)
                                .create();

    // Obtain kets from the two-atom bases and check for equality
    auto ket1a = combined_basis_a->get_kets()[0];
    auto ket1b = combined_basis_b->get_kets()[0];
    auto ket2a = combined_basis_a->get_kets()[1];
    auto ket2b = combined_basis_b->get_kets()[1];
    DOCTEST_CHECK(*ket1a == *ket1a);
    DOCTEST_CHECK(*ket2a == *ket2a);
    DOCTEST_CHECK(*ket1a != *ket2b);
    DOCTEST_CHECK(*ket2a != *ket1b);

    // Currently, kets from different combined bases are never equal
    DOCTEST_CHECK(*ket1a != *ket1b);
    DOCTEST_CHECK(*ket2a != *ket2b);
}
} // namespace pairinteraction
