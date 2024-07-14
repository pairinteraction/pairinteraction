#include "pintr/system/SystemAtom.hpp"

#include "pintr/basis/BasisAtom.hpp"
#include "pintr/basis/BasisAtomCreator.hpp"
#include "pintr/database/Database.hpp"
#include "pintr/diagonalizer/DiagonalizerFeast.hpp"
#include "pintr/diagonalizer/DiagonalizerLapacke.hpp"
#include "pintr/diagonalizer/diagonalize.hpp"
#include "pintr/ket/KetAtomCreator.hpp"

#include <Eigen/Eigenvalues>
#include <doctest/doctest.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

namespace pintr {
DOCTEST_TEST_CASE("construct and diagonalize a small Hamiltonian") {
    auto &database = Database::get_global_instance();
    auto diagonalizer = DiagonalizerLapacke<double>();

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
    auto diagonalizer = DiagonalizerLapacke<std::complex<double>>();

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
    auto diagonalizer = DiagonalizerLapacke<std::complex<double>>();

    auto basis = BasisAtomCreator<std::complex<double>>()
                     .set_species("Sr87_mqdt")
                     .restrict_quantum_number_nu(59, 61)
                     .restrict_quantum_number_l(0, 1)
                     .create(database);

    std::vector<SystemAtom<std::complex<double>>> systems;
    systems.reserve(n);
    for (int i = 0; i < n; ++i) {
        auto system = SystemAtom<std::complex<double>>(basis);
        system.set_electric_field({0, 0, 0.0001 * i});
        systems.emplace_back(system);
    }

    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Basis size: {}", basis->get_number_of_states());

    diagonalize<SystemAtom<std::complex<double>>>(systems, diagonalizer);
}

DOCTEST_TEST_CASE("construct and diagonalize a Hamiltonian using different methods") {
    auto &database = Database::get_global_instance();

    auto basis = BasisAtomCreator<std::complex<double>>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(59, 61)
                     .restrict_quantum_number_l(0, 1)
                     .create(database);

    for (int precision : {1, 2, 4, 8, 12}) {
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Precision: {}", precision);

        {
            auto diagonalizer = DiagonalizerLapacke<std::complex<double>>();

            auto system = SystemAtom<std::complex<double>>(basis);
            system.set_electric_field({0.0001, 0.0002, 0.0003});
            system.diagonalize(diagonalizer, precision);

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver;
            eigensolver.compute(system.get_matrix());
            auto eigenvalues_eigen = eigensolver.eigenvalues();
            auto eigenvalues_pairinteraction = system.get_matrix().diagonal();
            for (int i = 0; i < eigenvalues_eigen.size(); ++i) {
                DOCTEST_CHECK(std::abs(eigenvalues_eigen(i) - eigenvalues_pairinteraction(i)) <
                              10 * std::pow(10, -precision));
            }
        }

#ifdef WITH_MKL
        {
            auto diagonalizer = DiagonalizerFeast<std::complex<double>>(300);

            auto system = SystemAtom<std::complex<double>>(basis);
            system.set_electric_field({0.0001, 0.0002, 0.0003});
            system.diagonalize(diagonalizer, precision);

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver;
            eigensolver.compute(system.get_matrix());
            auto eigenvalues_eigen = eigensolver.eigenvalues();
            auto eigenvalues_pairinteraction = system.get_matrix().diagonal();
            for (int i = 0; i < eigenvalues_eigen.size(); ++i) {
                DOCTEST_CHECK(std::abs(eigenvalues_eigen(i) - eigenvalues_pairinteraction(i)) <
                              10 * std::pow(10, -precision));
            }
        }
#endif
    }
}
} // namespace pintr
