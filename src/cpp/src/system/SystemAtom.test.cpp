#include "pairinteraction/system/SystemAtom.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerEigen.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerFeast.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerLapacke.hpp"
#include "pairinteraction/diagonalizer/diagonalize.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"

#include <Eigen/Eigenvalues>
#include <doctest/doctest.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("construct and diagonalize a small Hamiltonian") {
    auto &database = Database::get_global_instance();
    auto diagonalizer = DiagonalizerEigen<double>();

    auto ket1 = KetAtomCreator<double>()
                    .set_species("Rb")
                    .set_quantum_number_n(60)
                    .set_quantum_number_l(0)
                    .set_quantum_number_j(0.5)
                    .set_quantum_number_m(0.5)
                    .create(database);
    auto ket2 = KetAtomCreator<double>()
                    .set_species("Rb")
                    .set_quantum_number_n(60)
                    .set_quantum_number_l(1)
                    .set_quantum_number_j(0.5)
                    .set_quantum_number_m(0.5)
                    .create(database);
    auto basis = BasisAtomCreator<double>().add_ket(ket1).add_ket(ket2).create(database);

    auto system = SystemAtom<double>(basis);
    system.set_electric_field({0, 0, 0.0001});

    Eigen::MatrixXd tmp = Eigen::MatrixXd(1e5 * system.get_matrix()).array().round() / 1e5;
    std::vector<double> matrix_vector(tmp.data(), tmp.data() + tmp.size());
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Constructed: {}", fmt::join(matrix_vector, ", "));

    system.diagonalize(diagonalizer);
    tmp = Eigen::MatrixXd(1e5 * system.get_matrix()).array().round() / 1e5;
    matrix_vector = std::vector<double>(tmp.data(), tmp.data() + tmp.size());
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Diagonalized: {}", fmt::join(matrix_vector, ", "));

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
    eigensolver.compute(system.get_matrix());
    auto eigenvalues_eigen = eigensolver.eigenvalues();
    auto eigenvalues_pairinteraction = system.get_matrix().diagonal();
    for (int i = 0; i < eigenvalues_eigen.size(); ++i) {
        DOCTEST_CHECK(std::abs(eigenvalues_eigen(i) - eigenvalues_pairinteraction(i)) < 1e-11);
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
                DOCTEST_CHECK(std::abs(matrix1.coeff(i, j)) < 1e-11);
                DOCTEST_CHECK(std::abs(matrix2.coeff(i, j)) < 1e-11);
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

    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Basis size: {}", basis->get_number_of_states());

    diagonalize<SystemAtom<std::complex<double>>>(systems, diagonalizer);
}

DOCTEST_TEST_CASE("construct and diagonalize a Hamiltonian using different methods") {
    auto &database = Database::get_global_instance();

    auto basis = BasisAtomCreator<std::complex<double>>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(60, 61)
                     .restrict_quantum_number_l(0, 1)
                     .create(database);

    std::vector<std::unique_ptr<DiagonalizerInterface<std::complex<double>>>> diagonalizers;
    diagonalizers.push_back(std::make_unique<DiagonalizerEigen<std::complex<double>>>());
#ifdef WITH_LAPACKE
    diagonalizers.push_back(std::make_unique<DiagonalizerLapacke<std::complex<double>>>());
#endif
#ifdef WITH_MKL
    diagonalizers.push_back(std::make_unique<DiagonalizerFeast<std::complex<double>>>(300));
#endif

    for (int precision : {1, 4, 12}) {
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Precision: {}", precision);

        for (const auto &diagonalizer : diagonalizers) {
            auto system = SystemAtom<std::complex<double>>(basis);
            system.set_electric_field({0.0001, 0.0002, 0.0003});

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver;
            eigensolver.compute(system.get_matrix());
            auto eigenvalues_eigen = eigensolver.eigenvalues();

            // We specify a search interval because this is required if the FEAST routine is used.
            // To avoid overflows, the interval ranges from half the smallest possible value to half
            // the largest possible value.
            system.diagonalize(*diagonalizer, 12,
                               {std::numeric_limits<double>::lowest() / 2,
                                std::numeric_limits<double>::max() / 2});
            auto eigenvalues_pairinteraction = system.get_matrix().diagonal();
            for (int i = 0; i < eigenvalues_eigen.size(); ++i) {
                DOCTEST_CHECK(std::abs(eigenvalues_eigen(i) - eigenvalues_pairinteraction(i)) <
                              10 * std::pow(10, -precision));
            }
        }
    }
}

DOCTEST_TEST_CASE("construct and diagonalize a Hamiltonian with energy restrictions") {
    double min_energy = -1.45e-4;
    double max_energy = -1.35e-4;

    auto &database = Database::get_global_instance();

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 1)
                     .create(database);

    std::vector<std::unique_ptr<DiagonalizerInterface<double>>> diagonalizers;
    diagonalizers.push_back(std::make_unique<DiagonalizerEigen<double>>());
#ifdef WITH_LAPACKE
    diagonalizers.push_back(std::make_unique<DiagonalizerLapacke<double>>());
#endif
#ifdef WITH_MKL
    diagonalizers.push_back(std::make_unique<DiagonalizerFeast<double>>(5));
#endif

    for (const auto &diagonalizer : diagonalizers) {
        auto system = SystemAtom<double>(basis);
        system.set_electric_field({0.0001, 0, 0.0001});

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
        eigensolver.compute(system.get_matrix());
        auto eigenvalues_all = eigensolver.eigenvalues();
        std::vector<double> eigenvalues_eigen;
        for (int i = 0; i < eigenvalues_all.size(); ++i) {
            if (eigenvalues_all[i] > min_energy && eigenvalues_all[i] < max_energy) {
                eigenvalues_eigen.push_back(eigenvalues_all[i]);
            }
        }

        system.diagonalize(*diagonalizer, 12, {min_energy, max_energy});
        auto eigenvalues_pairinteraction = system.get_matrix().diagonal();
        Eigen::MatrixXd tmp = (1e5 * eigenvalues_pairinteraction).array().round() / 1e5;
        std::vector<double> eigenvalues_vector(tmp.data(), tmp.data() + tmp.size());
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Eigenvalues: {}",
                           fmt::join(eigenvalues_vector, ", "));

        DOCTEST_CHECK(eigenvalues_eigen.size() == eigenvalues_pairinteraction.size());
        for (size_t i = 0; i < eigenvalues_eigen.size(); ++i) {
            DOCTEST_CHECK(std::abs(eigenvalues_eigen[i] - eigenvalues_pairinteraction[i]) < 1e-11);
        }
    }
}
} // namespace pairinteraction
