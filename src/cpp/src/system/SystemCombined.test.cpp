#include "pairinteraction/system/SystemCombined.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/basis/BasisCombined.hpp"
#include "pairinteraction/basis/BasisCombinedCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerEigen.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerFeast.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerLapacke.hpp"
#include "pairinteraction/diagonalizer/diagonalize.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/Range.hpp"

#include <Eigen/Eigenvalues>
#include <doctest/doctest.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("construct a combined Hamiltonian") {
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
    system1.set_electric_field({0, 0, 0.00000001});
    system2.set_electric_field({0, 0, 0.00000001});
    diagonalize<SystemAtom<double>>({system1, system2}, diagonalizer);

    // Construct and diagonalize the combined system
    auto basis_combined = BasisCombinedCreator<double>().add(system1).add(system2).create();

    auto system_combined = SystemCombined<double>(basis_combined);
    system_combined.set_distance(0.00000001);
    system_combined.diagonalize(diagonalizer);

    // Print the largest and smallest eigenvalues
    auto eigenvalues = system_combined.get_eigenvalues();
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Lowest energy: {}", eigenvalues.minCoeff());
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Highest energy: {}", eigenvalues.maxCoeff());
}
} // namespace pairinteraction
