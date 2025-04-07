// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/SystemPair.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/basis/BasisPairCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalize/DiagonalizerEigen.hpp"
#include "pairinteraction/diagonalize/DiagonalizerFeast.hpp"
#include "pairinteraction/diagonalize/DiagonalizerLapackeEvr.hpp"
#include "pairinteraction/diagonalize/diagonalize.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/Range.hpp"

#include <doctest/doctest.h>
#include <fmt/ranges.h>

namespace pairinteraction {

constexpr double VOLT_PER_CM_IN_ATOMIC_UNITS = 1 / 5.14220675112e9;
constexpr double UM_IN_ATOMIC_UNITS = 1 / 5.29177210544e-5;
constexpr double HARTREE_IN_GHZ = 6579683.920501762;

DOCTEST_TEST_CASE("construct a pair Hamiltonian") {
    auto &database = Database::get_global_instance();
    auto diagonalizer = DiagonalizerEigen<double>();

    // Construct and diagonalize the constituent systems
    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(60, 61)
                     .restrict_quantum_number_l(0, 1)
                     .restrict_quantum_number_m(-0.5, 0.5)
                     .create(database);

    SystemAtom<double> system1(basis);
    SystemAtom<double> system2(basis);
    system1.set_electric_field({0, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});
    system2.set_electric_field({0, 0, 2 * VOLT_PER_CM_IN_ATOMIC_UNITS});
    diagonalize<SystemAtom<double>>({system1, system2}, diagonalizer);

    // Construct and diagonalize the system_pair
    auto basis_pair = BasisPairCreator<double>().add(system1).add(system2).create();
    DOCTEST_MESSAGE("Number of states in pair basis: ", basis_pair->get_number_of_states());

    auto system_pair = SystemPair<double>(basis_pair);
    system_pair.set_distance_vector({0, 0, 3 * UM_IN_ATOMIC_UNITS});
    system_pair.diagonalize(diagonalizer);

    // Print the largest and smallest eigenenergies
    auto eigenenergies = system_pair.get_eigenenergies();
    DOCTEST_MESSAGE("Lowest energy: ", eigenenergies.minCoeff());
    DOCTEST_MESSAGE("Highest energy: ", eigenenergies.maxCoeff());
}

#ifdef WITH_LAPACKE
DOCTEST_TEST_CASE("diagonalize with lapacke_evr") {
    auto &database = Database::get_global_instance();
    auto diagonalizer = DiagonalizerLapackeEvr<double>();

    // Construct the state of interest
    auto ket = KetAtomCreator()
                   .set_species("Rb")
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(0)
                   .set_quantum_number_j(0.5)
                   .set_quantum_number_m(0.5)
                   .create(database);

    // Construct and diagonalize the constituent systems
    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(ket->get_quantum_number_n() - 3,
                                                ket->get_quantum_number_n() + 3)
                     .restrict_quantum_number_l(ket->get_quantum_number_l() - 1,
                                                ket->get_quantum_number_l() + 1)
                     .create(database);
    SystemAtom<double> system(basis);

    // Construct and diagonalize the system_pair
    auto basis_pair = BasisPairCreator<double>()
                          .add(system)
                          .add(system)
                          .restrict_energy(2 * ket->get_energy() - 2 / HARTREE_IN_GHZ,
                                           2 * ket->get_energy() + 2 / HARTREE_IN_GHZ)
                          .create();
    DOCTEST_MESSAGE("Number of states in pair basis: ", basis_pair->get_number_of_states());

    auto system_pair = SystemPair<double>(basis_pair);
    system_pair.set_distance_vector({0, 0, 3 * UM_IN_ATOMIC_UNITS});
    system_pair.diagonalize(diagonalizer, 2 * ket->get_energy() - 0.5 / HARTREE_IN_GHZ,
                            2 * ket->get_energy() + 0.5 / HARTREE_IN_GHZ);

    // Print the largest and smallest eigenenergies
    auto eigenenergies = system_pair.get_eigenenergies();
    DOCTEST_MESSAGE("Lowest energy: ", eigenenergies.minCoeff());
    DOCTEST_MESSAGE("Highest energy: ", eigenenergies.maxCoeff());
}
#endif
} // namespace pairinteraction
