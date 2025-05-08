// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
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
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"

#include <Eigen/Eigenvalues>
#include <doctest/doctest.h>
#include <fmt/ranges.h>

namespace pairinteraction {

constexpr double VOLT_PER_CM_IN_ATOMIC_UNITS = 1 / 5.14220675112e9;
constexpr double UM_IN_ATOMIC_UNITS = 1 / 5.29177210544e-5;
constexpr double HARTREE_IN_GHZ = 6579683.920501762;

DOCTEST_TEST_CASE("construct and diagonalize a small Hamiltonian") {
    auto &database = Database::get_global_instance();
    auto diagonalizer = DiagonalizerEigen<double>();

    auto ket1 = KetAtomCreator()
                    .set_species("Rb")
                    .set_quantum_number_n(60)
                    .set_quantum_number_l(0)
                    .set_quantum_number_j(0.5)
                    .set_quantum_number_m(0.5)
                    .create(database);
    auto ket2 = KetAtomCreator()
                    .set_species("Rb")
                    .set_quantum_number_n(60)
                    .set_quantum_number_l(1)
                    .set_quantum_number_j(0.5)
                    .set_quantum_number_m(0.5)
                    .create(database);
    auto basis = BasisAtomCreator<double>().append_ket(ket1).append_ket(ket2).create(database);

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
                     .restrict_quantum_number_n(59, 61)
                     .restrict_quantum_number_l(0, 1)
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
                     .restrict_quantum_number_nu(60, 61)
                     .restrict_quantum_number_l(0, 1)
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
                     .restrict_quantum_number_n(60, 61)
                     .restrict_quantum_number_l(0, 1)
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
            std::max(0.5 * rtol_eigenenergies / std::sqrt(basis->get_number_of_states()), 4 * eps);
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
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 1)
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
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 1)
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
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(1)
                   .set_quantum_number_j(0.5)
                   .set_quantum_number_m(0.5)
                   .create(database);
    double energy = ket->get_energy();
    double min_energy = energy - 50 / HARTREE_IN_GHZ;
    double max_energy = energy + 50 / HARTREE_IN_GHZ;

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 3)
                     .restrict_quantum_number_m(0.5, 0.5)
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

} // namespace pairinteraction
