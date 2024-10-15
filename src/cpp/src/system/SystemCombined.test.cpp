#include "pairinteraction/system/SystemCombined.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
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

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(60, 61)
                     .restrict_quantum_number_l(0, 1)
                     .restrict_quantum_number_m(-0.5, 0.5)
                     .create(database);

    // Construct and diagonalize the constituent systems
    auto system1 = SystemAtom<double>(basis);
    auto system2 = SystemAtom<double>(basis);
    system1.set_electric_field({0, 0, 0.00000001});
    system2.set_electric_field({0, 0, 0.00000001});
    diagonalize<SystemAtom<double>>({system1, system2}, diagonalizer);

    // Construct and diagonalize the combined system
    auto combined = SystemCombined<double>(system1, system2, -1, 0);
    combined.set_interatomic_distance(0.00000001);
    combined.diagonalize(diagonalizer);

    // Print the largest and smallest eigenvalues
    auto eigenvalues = combined.get_eigenvalues();
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Lowest energy: {}", eigenvalues.minCoeff());
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Highest energy: {}", eigenvalues.maxCoeff());
}
} // namespace pairinteraction
