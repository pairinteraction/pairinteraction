// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/SystemAtom.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalize/DiagonalizerEigen.hpp"
#include "pairinteraction/diagonalize/DiagonalizerFeast.hpp"
#include "pairinteraction/diagonalize/DiagonalizerLapackeEvd.hpp"
#include "pairinteraction/diagonalize/DiagonalizerLapackeEvr.hpp"
#include "pairinteraction/diagonalize/diagonalize.hpp"
#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/SorterType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <doctest/doctest.h>
#include <fmt/ranges.h>
#include <stdexcept>

namespace pairinteraction {

constexpr double VOLT_PER_CM_IN_ATOMIC_UNITS = 1 / 5.14220675112e9;
constexpr double UM_IN_ATOMIC_UNITS = 1 / 5.29177210544e-5;
constexpr double HARTREE_IN_GHZ = 6579683.920501762;

DOCTEST_TEST_CASE("construct and diagonalize a small Hamiltonian") {
    auto &database = Database::get_global_instance();
    auto diagonalizer = DiagonalizerEigen<double>();

    auto ket1 = KetAtomCreator()
                    .set_species("Rb")
                    .set_quantum_number("n", 60)
                    .set_quantum_number("l", 0)
                    .set_quantum_number("j", 0.5)
                    .set_quantum_number("m", 0.5)
                    .create(database);
    auto ket2 = KetAtomCreator()
                    .set_species("Rb")
                    .set_quantum_number("n", 60)
                    .set_quantum_number("l", 1)
                    .set_quantum_number("j", 0.5)
                    .set_quantum_number("m", 0.5)
                    .create(database);
    auto basis = BasisAtomCreator<double>().add_ket(ket1).add_ket(ket2).create(database);

    auto system = SystemAtom<double>(basis);
    system.set_electric_field({0, 0, 0.0001});

    Eigen::MatrixXd tmp = Eigen::MatrixXd(1e5 * system.get_matrix()).array().round() / 1e5;
    std::vector<double> matrix_vector(tmp.data(), tmp.data() + tmp.size());
    DOCTEST_MESSAGE(fmt::format("Constructed: {}", fmt::join(matrix_vector, ", ")));

    system.diagonalize(diagonalizer);
    tmp = Eigen::MatrixXd(1e5 * system.get_matrix()).array().round() / 1e5;
    matrix_vector = std::vector<double>(tmp.data(), tmp.data() + tmp.size());
    DOCTEST_MESSAGE(fmt::format("Diagonalized: {}", fmt::join(matrix_vector, ", ")));

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
    eigensolver.compute(system.get_matrix());
    auto eigenenergies_eigen = eigensolver.eigenvalues();
    auto eigenenergies_pairinteraction = system.get_eigenenergies();
    for (int i = 0; i < eigenenergies_eigen.size(); ++i) {
        DOCTEST_CHECK(std::abs(eigenenergies_eigen(i) - eigenenergies_pairinteraction(i)) < 1e-10);
    }
}

DOCTEST_TEST_CASE("construct and diagonalize two Hamiltonians in parallel") {
    auto &database = Database::get_global_instance();
    auto diagonalizer = DiagonalizerEigen<std::complex<double>>();

    auto basis = BasisAtomCreator<std::complex<double>>()
                     .set_species("Rb")
                     .restrict_quantum_number("n", 59, 61)
                     .restrict_quantum_number("l", 0, 1)
                     .create(database);

    auto system1 = SystemAtom<std::complex<double>>(basis);
    system1.set_electric_field({0, 0, 0.0001});

    auto system2 = SystemAtom<std::complex<double>>(basis);
    system2.set_electric_field({0, 0, 0.0002});

    diagonalize<SystemAtom<std::complex<double>>>({system1, system2}, diagonalizer);

    auto matrix1 = system1.get_matrix();
    auto matrix2 = system2.get_matrix();
    for (int i = 0; i < matrix1.rows(); ++i) {
        for (int j = 0; j < matrix1.cols(); ++j) {
            if (i != j) {
                DOCTEST_CHECK(std::abs(matrix1.coeff(i, j)) < 1e-10);
                DOCTEST_CHECK(std::abs(matrix2.coeff(i, j)) < 1e-10);
            }
        }
    }
}

DOCTEST_TEST_CASE("construct an atomic Hamiltonian in a non-canonical atomic basis") {
    auto &database = Database::get_global_instance();

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number("n", 60, 61)
                     .restrict_quantum_number("l", 0, 1)
                     .restrict_quantum_number("m", 0.5, 0.5)
                     .create(database);
    DOCTEST_REQUIRE(basis->get_number_of_states() >= 2);

    SystemAtom<double> reference_system(basis);
    reference_system.set_electric_field({0, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});
    const auto &reference_matrix = reference_system.get_matrix();

    Eigen::SparseMatrix<double, Eigen::RowMajor> transformation(
        static_cast<Eigen::Index>(basis->get_number_of_states()),
        static_cast<Eigen::Index>(basis->get_number_of_states()));
    transformation.setIdentity();

    double inverse_sqrt_two = 1 / std::sqrt(2.0);
    transformation.coeffRef(0, 0) = inverse_sqrt_two;
    transformation.coeffRef(1, 0) = inverse_sqrt_two;
    transformation.coeffRef(0, 1) = inverse_sqrt_two;
    transformation.coeffRef(1, 1) = -inverse_sqrt_two;
    transformation.makeCompressed();

    auto transformed_basis = basis->transformed(transformation);
    SystemAtom<double> transformed_system(transformed_basis);
    transformed_system.set_electric_field({0, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});

    Eigen::SparseMatrix<double, Eigen::RowMajor> expected_matrix =
        transformation.adjoint() * reference_matrix * transformation;

    DOCTEST_CHECK(transformed_system.get_matrix().isApprox(expected_matrix, 1e-11));
}

DOCTEST_TEST_CASE("construct and diagonalize multiple Hamiltonians in parallel" *
                  doctest::skip(true)) {
    // TODO For a slow database, the fast parallelized construction of the tiny Hamiltonians seems
    // to lead to the wrong Hamiltonians (visible in the warning "The floating point error (5e-324
    // Hartree) is similar or larger than error estimated from the specified tolerance (0
    // Hartree)."). This could be caused by an issue with thread safety or memory access.
    int n = 10;

    auto &database = Database::get_global_instance();
    auto diagonalizer = DiagonalizerEigen<std::complex<double>>();

    auto basis = BasisAtomCreator<std::complex<double>>()
                     .set_species("Sr87_mqdt")
                     .restrict_quantum_number("nu", 60, 61)
                     .restrict_quantum_number("l", 0, 1)
                     .create(database);

    std::vector<SystemAtom<std::complex<double>>> systems;
    systems.reserve(n);
    for (int i = 0; i < n; ++i) {
        auto system = SystemAtom<std::complex<double>>(basis);
        system.set_electric_field({0, 0, 0.0001 * i});
        systems.push_back(std::move(system));
    }

    DOCTEST_MESSAGE("Basis size: ", basis->get_number_of_states());

    diagonalize<SystemAtom<std::complex<double>>>(systems, diagonalizer);
}

DOCTEST_TEST_CASE("construct and diagonalize a Hamiltonian using different methods") {
    auto &database = Database::get_global_instance();

    auto basis = BasisAtomCreator<std::complex<double>>()
                     .set_species("Rb")
                     .restrict_quantum_number("n", 60, 61)
                     .restrict_quantum_number("l", 0, 1)
                     .create(database);

    // Diagonalize using the Eigen library
    auto system = SystemAtom<std::complex<double>>(basis);
    system.set_electric_field({1 * VOLT_PER_CM_IN_ATOMIC_UNITS, 2 * VOLT_PER_CM_IN_ATOMIC_UNITS,
                               3 * VOLT_PER_CM_IN_ATOMIC_UNITS});

    Eigen::MatrixXcd matrix = system.get_matrix();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver;
    eigensolver.compute(matrix);
    auto eigenenergies_eigen = eigensolver.eigenvalues();

    // Create diagonalizers
    std::vector<std::unique_ptr<DiagonalizerInterface<std::complex<double>>>> diagonalizers;
    std::vector<double> rtols;
    double eps{};
    DOCTEST_SUBCASE("Double precision") {
        diagonalizers.push_back(std::make_unique<DiagonalizerEigen<std::complex<double>>>());
#ifdef WITH_LAPACKE
        diagonalizers.push_back(std::make_unique<DiagonalizerLapackeEvd<std::complex<double>>>());
        diagonalizers.push_back(std::make_unique<DiagonalizerLapackeEvr<std::complex<double>>>());
#endif
#ifdef WITH_MKL
        diagonalizers.push_back(std::make_unique<DiagonalizerFeast<std::complex<double>>>(300));
#endif
        rtols = {1e-1, 1e-6, 1e-14};
        eps = std::numeric_limits<double>::epsilon();
    }

    DOCTEST_SUBCASE("Single precision") {
        diagonalizers.push_back(
            std::make_unique<DiagonalizerEigen<std::complex<double>>>(FloatType::FLOAT32));
#ifdef WITH_LAPACKE
        diagonalizers.push_back(
            std::make_unique<DiagonalizerLapackeEvd<std::complex<double>>>(FloatType::FLOAT32));
        diagonalizers.push_back(
            std::make_unique<DiagonalizerLapackeEvr<std::complex<double>>>(FloatType::FLOAT32));
#endif
#ifdef WITH_MKL
        diagonalizers.push_back(
            std::make_unique<DiagonalizerFeast<std::complex<double>>>(300, FloatType::FLOAT32));
#endif
        rtols = {1e-1, 1e-6};
        eps = std::numeric_limits<float>::epsilon();
    }

    // Diagonalize using pairinteraction
    for (double rtol_eigenenergies : rtols) {
        double atol_eigenvectors =
            std::max(0.5 * rtol_eigenenergies / std::sqrt(basis->get_number_of_states()), 5 * eps);
        DOCTEST_MESSAGE("Precision: " << rtol_eigenenergies << " (rtol eigenenergies), "
                                      << atol_eigenvectors << " (atol eigenvectors)");

        for (const auto &diagonalizer : diagonalizers) {
            auto system = SystemAtom<std::complex<double>>(basis);
            system.set_electric_field({1 * VOLT_PER_CM_IN_ATOMIC_UNITS,
                                       2 * VOLT_PER_CM_IN_ATOMIC_UNITS,
                                       3 * VOLT_PER_CM_IN_ATOMIC_UNITS});

            // We specify a search interval because this is required if the FEAST routine is
            // used. To avoid overflows, the interval ranges from half the smallest possible
            // value to half the largest possible value.
            system.diagonalize(*diagonalizer, std::numeric_limits<float>::lowest() / 2,
                               std::numeric_limits<float>::max() / 2, rtol_eigenenergies);
            auto eigenenergies_pairinteraction = system.get_eigenenergies();
            auto eigenvectors_pairinteraction = system.get_eigenbasis()->get_coefficients();

            DOCTEST_CHECK(
                (eigenenergies_eigen - eigenenergies_pairinteraction).array().abs().maxCoeff() <
                rtol_eigenenergies * matrix.norm());

            for (int i = 0; i < eigenvectors_pairinteraction.cols(); ++i) {
                DOCTEST_CHECK(abs(1 - eigenvectors_pairinteraction.col(i).norm()) <
                              atol_eigenvectors * eigenvectors_pairinteraction.rows());
            }
        }
    }
}

DOCTEST_TEST_CASE("construct and diagonalize a Hamiltonian with energy restrictions") {
    double min_energy = 0.153355;
    double max_energy = 0.153360;

    auto &database = Database::get_global_instance();

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number("n", 58, 62)
                     .restrict_quantum_number("l", 0, 1)
                     .create(database);

    // Diagonalize using the Eigen library
    auto system = SystemAtom<double>(basis);
    system.set_electric_field(
        {1 * VOLT_PER_CM_IN_ATOMIC_UNITS, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
    eigensolver.compute(system.get_matrix());
    auto eigenenergies_all = eigensolver.eigenvalues();
    std::vector<double> eigenenergies_eigen;
    for (int i = 0; i < eigenenergies_all.size(); ++i) {
        if (eigenenergies_all[i] > min_energy && eigenenergies_all[i] < max_energy) {
            eigenenergies_eigen.push_back(eigenenergies_all[i]);
        }
    }

    // Create diagonalizer
    std::vector<std::unique_ptr<DiagonalizerInterface<double>>> diagonalizers;
    diagonalizers.push_back(std::make_unique<DiagonalizerEigen<double>>(FloatType::FLOAT64));
#ifdef WITH_LAPACKE
    diagonalizers.push_back(std::make_unique<DiagonalizerLapackeEvd<double>>(FloatType::FLOAT64));
    diagonalizers.push_back(std::make_unique<DiagonalizerLapackeEvr<double>>(FloatType::FLOAT64));
#endif
#ifdef WITH_MKL
    diagonalizers.push_back(std::make_unique<DiagonalizerFeast<double>>(10, FloatType::FLOAT64));
#endif

    // Diagonalize using pairinteraction
    for (const auto &diagonalizer : diagonalizers) {
        auto system = SystemAtom<double>(basis);
        system.set_electric_field(
            {1 * VOLT_PER_CM_IN_ATOMIC_UNITS, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});

        system.diagonalize(*diagonalizer, min_energy, max_energy, 1e-6);
        auto eigenenergies_pairinteraction = system.get_eigenenergies();

        Eigen::MatrixXd tmp = (1e5 * eigenenergies_pairinteraction).array().round() / 1e5;
        std::vector<double> eigenenergies_vector(tmp.data(), tmp.data() + tmp.size());
        DOCTEST_MESSAGE(fmt::format("Eigenenergies: {}", fmt::join(eigenenergies_vector, ", ")));

        DOCTEST_CHECK(eigenenergies_eigen.size() == 8);
        DOCTEST_CHECK(eigenenergies_pairinteraction.size() == 8);
        for (size_t i = 0; i < eigenenergies_eigen.size(); ++i) {
            DOCTEST_CHECK(std::abs(eigenenergies_eigen[i] - eigenenergies_pairinteraction[i]) <
                          1e-10);
        }
    }
}

#ifdef WITH_MKL
#include <Eigen/Dense>
#include <mkl.h>
DOCTEST_TEST_CASE("diagonalization with mkl") {
    // We loop several times to check for memory errors
    for (size_t i = 0; i < 10; ++i) {
        // Create a symmetric matrix
        int n = 100;
        Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(n, n);
        matrix = (matrix + matrix.transpose()).eval();
        Eigen::VectorXd eigenenergies(n);

        // Diagonalize the matrix
        int info =
            LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, matrix.data(), n, eigenenergies.data());
        DOCTEST_CHECK(info == 0);
    }
}
#endif

DOCTEST_TEST_CASE("handle it gracefully if no eigenenergies are within energy restrictions") {
    double min_energy = -1;
    double max_energy = -1;

    auto &database = Database::get_global_instance();

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number("n", 58, 62)
                     .restrict_quantum_number("l", 0, 1)
                     .create(database);

    std::vector<std::unique_ptr<DiagonalizerInterface<double>>> diagonalizers;
    diagonalizers.push_back(std::make_unique<DiagonalizerEigen<double>>());
#ifdef WITH_LAPACKE
    diagonalizers.push_back(std::make_unique<DiagonalizerLapackeEvd<double>>());
#endif

    for (const auto &diagonalizer : diagonalizers) {
        auto system = SystemAtom<double>(basis);
        system.set_electric_field(
            {1 * VOLT_PER_CM_IN_ATOMIC_UNITS, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});

        system.diagonalize(*diagonalizer, min_energy, max_energy, 1e-6);
        auto eigenenergies_pairinteraction = system.get_eigenenergies();

        DOCTEST_CHECK(eigenenergies_pairinteraction.size() == 0);
    }
}

DOCTEST_TEST_CASE("atom ion interaction") {
    auto &database = Database::get_global_instance();
    DiagonalizerEigen<double> diagonalizer;

    auto ket = KetAtomCreator()
                   .set_species("Rb")
                   .set_quantum_number("n", 60)
                   .set_quantum_number("l", 1)
                   .set_quantum_number("j", 0.5)
                   .set_quantum_number("m", 0.5)
                   .create(database);
    double energy = ket->get_energy();
    double min_energy = energy - 50 / HARTREE_IN_GHZ;
    double max_energy = energy + 50 / HARTREE_IN_GHZ;

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number("n", 58, 62)
                     .restrict_quantum_number("l", 0, 3)
                     .restrict_quantum_number("m", 0.5, 0.5)
                     .create(database);

    auto system3 = SystemAtom<double>(basis);
    system3.set_ion_interaction_order(3);
    system3.set_ion_distance_vector({0, 0, 10 * UM_IN_ATOMIC_UNITS});
    system3.diagonalize(diagonalizer, min_energy, max_energy, 1e-6);
    auto energies3 = system3.get_eigenenergies();

    auto system2 = SystemAtom<double>(basis);
    system2.set_ion_interaction_order(2);
    system2.set_ion_distance_vector({0, 0, 10 * UM_IN_ATOMIC_UNITS});
    system2.diagonalize(diagonalizer, min_energy, max_energy, 1e-6);
    auto energies2 = system2.get_eigenenergies();

    // Ensure that the quadrupole order has a significant effect
    size_t num_energies = std::min(energies2.size(), energies3.size());
    for (size_t i = 0; i < num_energies; ++i) {
        DOCTEST_CHECK(std::abs(energies3[i] - energies2[i]) * HARTREE_IN_GHZ > 1e-6);
    }
}

DOCTEST_TEST_CASE("obtain the blocks of a Hamiltonian") {
    auto &database = Database::get_global_instance();

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number("n", 60, 60)
                     .restrict_quantum_number("l", 0, 1)
                     .create(database);

    auto system = SystemAtom<double>(basis);
    system.set_electric_field({0, 0, 0.0001});

    // The electric field conserves the quantum number m, thus the Hamiltonian can be
    // block-diagonalized by m
    system.transform(system.get_sorter({SorterType::QUANTUM_NUMBER_M}));
    auto blocks = system.get_indices_of_blocks({SorterType::QUANTUM_NUMBER_M});
    DOCTEST_CHECK(blocks.size() > 1);

    size_t expected_start = 0;
    for (const auto &block : blocks) {
        DOCTEST_MESSAGE("Block from ", block.start, " to ", block.end);
        DOCTEST_CHECK(block.start == expected_start);
        expected_start = block.end;
    }
    DOCTEST_CHECK(expected_start == basis->get_number_of_states());

    // Blocks are obtained to diagonalize the Hamiltonian block by block. Because energy blocks
    // would require an already diagonal Hamiltonian, they cannot be obtained ...
    DOCTEST_CHECK_THROWS_AS(system.get_indices_of_blocks({SorterType::ENERGY}),
                            std::invalid_argument);

    // ... whereas sorting by the energy is supported
    DOCTEST_CHECK_NOTHROW(system.get_sorter({SorterType::ENERGY}));

    // A repeated label cannot influence the order anymore and must be ignored instead of
    // terminating the program
    DOCTEST_CHECK_NOTHROW(system.get_sorter({SorterType::ENERGY, SorterType::ENERGY}));
    DOCTEST_CHECK_NOTHROW(
        system.get_sorter({SorterType::QUANTUM_NUMBER_M, SorterType::ENERGY, SorterType::ENERGY}));
}

DOCTEST_TEST_CASE("sort a Hamiltonian by several labels") {
    auto &database = Database::get_global_instance();

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number("n", 60, 60)
                     .restrict_quantum_number("l", 0, 1)
                     .create(database);

    // The first label is the primary sorting criterion, the following labels only break ties
    {
        auto system = SystemAtom<double>(basis);
        system.set_electric_field({0, 0, 0.0001});
        system.transform(system.get_sorter({SorterType::QUANTUM_NUMBER_M, SorterType::ENERGY}));

        // Sorting by m first means that states of equal m are contiguous ...
        auto blocks = system.get_indices_of_blocks({SorterType::QUANTUM_NUMBER_M});
        DOCTEST_CHECK(blocks.size() > 1);

        // ... and that within each block, the energies are ascending
        const auto &matrix = system.get_matrix();
        for (const auto &block : blocks) {
            for (auto i = static_cast<long>(block.start) + 1; i < static_cast<long>(block.end);
                 ++i) {
                DOCTEST_CHECK(matrix.coeff(i - 1, i - 1) <= matrix.coeff(i, i));
            }
        }
    }

    // If the energy is the first label, the states are sorted by the energy globally
    {
        auto system = SystemAtom<double>(basis);
        system.set_electric_field({0, 0, 0.0001});
        system.transform(system.get_sorter({SorterType::ENERGY, SorterType::QUANTUM_NUMBER_M}));

        const auto &matrix = system.get_matrix();
        for (long i = 1; i < matrix.rows(); ++i) {
            DOCTEST_CHECK(matrix.coeff(i - 1, i - 1) <= matrix.coeff(i, i));
        }
    }

    // If the energy is enclosed by other labels, the labels before the energy take precedence over
    // the energy and the energy takes precedence over the labels after it
    {
        auto system = SystemAtom<double>(basis);
        system.set_electric_field({0, 0, 0.0001});
        system.transform(system.get_sorter(
            {SorterType::PARITY, SorterType::ENERGY, SorterType::QUANTUM_NUMBER_M}));

        auto sorted_basis = system.get_basis();
        size_t num_states = sorted_basis->get_number_of_states();
        const auto &matrix = system.get_matrix();

        auto parity = [&](size_t i) { return sorted_basis->get_quantum_number("parity", i); };
        auto energy = [&](size_t i) {
            return matrix.coeff(static_cast<long>(i), static_cast<long>(i));
        };
        auto quantum_number_m = [&](size_t i) { return sorted_basis->get_quantum_number("m", i); };

        // The parity is the primary criterion, thus it is ascending globally ...
        for (size_t i = 1; i < num_states; ++i) {
            DOCTEST_CHECK(parity(i - 1) <= parity(i));
        }

        // ... the energy is only ascending within a block of equal parity ...
        for (size_t i = 1; i < num_states; ++i) {
            if (parity(i - 1) == parity(i)) {
                DOCTEST_CHECK(energy(i - 1) <= energy(i));
            }
        }

        // ... and the quantum number m only breaks ties between states of equal parity and energy
        size_t num_tie_breaks = 0;
        for (size_t i = 1; i < num_states; ++i) {
            if (parity(i - 1) == parity(i) && energy(i - 1) == energy(i)) {
                DOCTEST_CHECK(quantum_number_m(i - 1) <= quantum_number_m(i));
                ++num_tie_breaks;
            }
        }

        // Check that the basis contains degenerate states so that the tie breaking by m above is
        // not vacuous
        DOCTEST_CHECK(num_tie_breaks > 0);

        // If m was the primary criterion, m would be ascending globally instead
        bool m_ascending_globally = true;
        for (size_t i = 1; i < num_states; ++i) {
            if (quantum_number_m(i - 1) > quantum_number_m(i)) {
                m_ascending_globally = false;
            }
        }
        DOCTEST_CHECK_FALSE(m_ascending_globally);
    }
}

} // namespace pairinteraction
