#include "pairinteraction/system/SystemPair.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/basis/BasisPairCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalize/DiagonalizerEigen.hpp"
#include "pairinteraction/diagonalize/DiagonalizerFeast.hpp"
#include "pairinteraction/diagonalize/DiagonalizerLapackeEvd.hpp"
#include "pairinteraction/diagonalize/diagonalize.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/Range.hpp"

#include <Eigen/Eigenvalues>
#include <doctest/doctest.h>
#include <fmt/ranges.h>

namespace pairinteraction {

constexpr double VOLT_PER_CM_IN_ATOMIC_UNITS = 1 / 5.14220675112e9;
constexpr double UM_IN_ATOMIC_UNITS = 1 / 5.29177210544e-5;

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

    auto system_pair = SystemPair<double>(basis_pair);
    system_pair.set_distance_vector({0, 0, 3 * UM_IN_ATOMIC_UNITS});
    system_pair.diagonalize(diagonalizer);

    // Print the largest and smallest eigenvalues
    auto eigenvalues = system_pair.get_eigenvalues();
    DOCTEST_MESSAGE("Lowest energy: ", eigenvalues.minCoeff());
    DOCTEST_MESSAGE("Highest energy: ", eigenvalues.maxCoeff());
}
} // namespace pairinteraction
