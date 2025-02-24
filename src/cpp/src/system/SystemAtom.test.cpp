#include "pairinteraction/system/SystemAtom.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerEigen.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerFeast.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerLapackeEvd.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerLapackeEvr.hpp"
#include "pairinteraction/diagonalizer/diagonalize.hpp"
#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"

#include <Eigen/Eigenvalues>
#include <doctest/doctest.h>
#include <fmt/ranges.h>

namespace pairinteraction {

constexpr double VOLT_PER_CM_IN_ATOMIC_UNITS = 1 / 5.14220675112e9;

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
    auto eigenvalues_eigen = eigensolver.eigenvalues();
    auto eigenvalues_pairinteraction = system.get_eigenvalues();
    for (int i = 0; i < eigenvalues_eigen.size(); ++i) {
        DOCTEST_CHECK(std::abs(eigenvalues_eigen(i) - eigenvalues_pairinteraction(i)) < 1e-10);
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

DOCTEST_TEST_CASE("construct and diagonalize multiple Hamiltonians in parallel") {
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
    auto eigenvalues_eigen = eigensolver.eigenvalues();
    auto eigenvalues_eigen_shifted = eigenvalues_eigen;
    eigenvalues_eigen_shifted.array() -= matrix.diagonal().real().mean();

    // Create diagonalizers
    std::vector<std::unique_ptr<DiagonalizerInterface<std::complex<double>>>> diagonalizers;
    std::vector<double> atols;
    DOCTEST_SUBCASE("Double precision") {
        diagonalizers.push_back(std::make_unique<DiagonalizerEigen<std::complex<double>>>());
#ifdef WITH_LAPACKE
        diagonalizers.push_back(std::make_unique<DiagonalizerLapackeEvd<std::complex<double>>>());
        diagonalizers.push_back(std::make_unique<DiagonalizerLapackeEvr<std::complex<double>>>());
#endif
#ifdef WITH_MKL
        diagonalizers.push_back(std::make_unique<DiagonalizerFeast<std::complex<double>>>(300));
#endif
        atols = {1e-1, 1e-6, 1e-14};
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
        atols = {1e-1, 1e-6};
    }

    // Diagonalize using pairinteraction
    for (double atol_eigenvectors : atols) {
        double atol_eigenvalues =
            std::max(1e-16, 2 * atol_eigenvectors * eigenvalues_eigen_shifted.norm());
        DOCTEST_MESSAGE("Precision: " << atol_eigenvectors << " (eigenvectors), "
                                      << atol_eigenvalues << " Hartree (eigenvalues)");

        for (const auto &diagonalizer : diagonalizers) {
            auto system = SystemAtom<std::complex<double>>(basis);
            system.set_electric_field({1 * VOLT_PER_CM_IN_ATOMIC_UNITS,
                                       2 * VOLT_PER_CM_IN_ATOMIC_UNITS,
                                       3 * VOLT_PER_CM_IN_ATOMIC_UNITS});

            // We specify a search interval because this is required if the FEAST routine is
            // used. To avoid overflows, the interval ranges from half the smallest possible
            // value to half the largest possible value.
            system.diagonalize(*diagonalizer, std::numeric_limits<float>::lowest() / 2,
                               std::numeric_limits<float>::max() / 2, atol_eigenvectors);
            auto eigenvalues_pairinteraction = system.get_eigenvalues();
            auto eigenvectors_pairinteraction = system.get_eigenbasis()->get_coefficients();

            DOCTEST_CHECK((eigenvalues_eigen - eigenvalues_pairinteraction).array().abs().sum() <
                          2 * atol_eigenvalues * eigenvalues_eigen.size()); // 2 is a safety factor

            for (int i = 0; i < eigenvectors_pairinteraction.cols(); ++i) {
                DOCTEST_CHECK(abs(1 - eigenvectors_pairinteraction.col(i).norm()) < 2 *
                                  atol_eigenvectors *
                                  eigenvectors_pairinteraction.rows()); // 2 is a safety factor
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
    auto eigenvalues_all = eigensolver.eigenvalues();
    std::vector<double> eigenvalues_eigen;
    for (int i = 0; i < eigenvalues_all.size(); ++i) {
        if (eigenvalues_all[i] > min_energy && eigenvalues_all[i] < max_energy) {
            eigenvalues_eigen.push_back(eigenvalues_all[i]);
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
        auto eigenvalues_pairinteraction = system.get_eigenvalues();

        Eigen::MatrixXd tmp = (1e5 * eigenvalues_pairinteraction).array().round() / 1e5;
        std::vector<double> eigenvalues_vector(tmp.data(), tmp.data() + tmp.size());
        DOCTEST_MESSAGE(fmt::format("Eigenvalues: {}", fmt::join(eigenvalues_vector, ", ")));

        DOCTEST_CHECK(eigenvalues_eigen.size() == 8);
        DOCTEST_CHECK(eigenvalues_pairinteraction.size() == 8);
        for (size_t i = 0; i < eigenvalues_eigen.size(); ++i) {
            DOCTEST_CHECK(std::abs(eigenvalues_eigen[i] - eigenvalues_pairinteraction[i]) < 1e-10);
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
        Eigen::VectorXd eigenvalues(n);

        // Diagonalize the matrix
        int info =
            LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, matrix.data(), n, eigenvalues.data());
        DOCTEST_CHECK(info == 0);
    }
}
#endif

DOCTEST_TEST_CASE("handle it gracefully if no eigenvalues are within energy restrictions") {
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
        auto eigenvalues_pairinteraction = system.get_eigenvalues();

        DOCTEST_CHECK(eigenvalues_pairinteraction.size() == 0);
    }
}

} // namespace pairinteraction
