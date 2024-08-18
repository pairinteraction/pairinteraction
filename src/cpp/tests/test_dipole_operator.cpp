#include "pairinteraction/pairinteraction.hpp"
#include "pairinteraction/utils/args.hpp"

#include <filesystem>
#include <spdlog/spdlog.h>

int main(int argc, char **argv) {
    // Call the setup function to configure logging
    pairinteraction::setup();

    // Create a database instance
    std::filesystem::path databasedir;
    bool download_missing = false;

    for (int i = 1; i < argc; ++i) {
        bool found = pairinteraction::args::parse_download_missing(i, argc, argv, download_missing);
        if (!found) {
            pairinteraction::args::parse_database(i, argc, argv, databasedir);
        }
    }

    pairinteraction::Database database(download_missing, true, databasedir);

    // Create a dipole operator coupling two specific states
    auto ket1 = pairinteraction::KetAtomCreator<float>()
                    .set_species("Rb")
                    .set_quantum_number_n(60)
                    .set_quantum_number_l(0)
                    .set_quantum_number_j(0.5)
                    .set_quantum_number_m(0.5)
                    .create(database);

    auto ket2 = pairinteraction::KetAtomCreator<float>()
                    .set_species("Rb")
                    .set_quantum_number_n(60)
                    .set_quantum_number_l(1)
                    .set_quantum_number_j(0.5)
                    .set_quantum_number_m(0.5)
                    .create(database);

    auto basis_ket1_ket2 =
        pairinteraction::BasisAtomCreator<float>().add_ket(ket1).add_ket(ket2).create(database);

    pairinteraction::OperatorAtom<float> dipole_ket1_ket2(
        basis_ket1_ket2, pairinteraction::OperatorType::ELECTRIC_DIPOLE, 0);
    float dipole_ket1_ket2_value = dipole_ket1_ket2.get_matrix().coeff(0, 1);

    if (std::abs(dipole_ket1_ket2_value - 1247.5955810546875) >
        10 * std::numeric_limits<float>::epsilon()) {
        SPDLOG_ERROR("The dipole operator value is not correct.");
        return 1;
    }

    dipole_ket1_ket2 = 2 * dipole_ket1_ket2;
    dipole_ket1_ket2_value = dipole_ket1_ket2.get_matrix().coeff(0, 1);

    if (std::abs(dipole_ket1_ket2_value - 2 * 1247.5955810546875) >
        10 * std::numeric_limits<float>::epsilon()) {
        SPDLOG_ERROR("The dipole operator value is not correct after multiplication by a scalar.");
        return 1;
    }

    dipole_ket1_ket2 = dipole_ket1_ket2 + dipole_ket1_ket2;
    dipole_ket1_ket2_value = dipole_ket1_ket2.get_matrix().coeff(0, 1);

    if (std::abs(dipole_ket1_ket2_value - 4 * 1247.5955810546875) >
        10 * std::numeric_limits<float>::epsilon()) {
        SPDLOG_ERROR("The dipole operator value is not correct after summation.");
        return 1;
    }

    // Create dipole operators in a typical basis
    auto basis = pairinteraction::BasisAtomCreator<std::complex<float>>()
                     .set_species("Sr88_singlet")
                     .restrict_quantum_number_n(60, 63)
                     .restrict_quantum_number_l(0, 3)
                     .create(database);

    pairinteraction::OperatorAtom<std::complex<float>> dipole_0(
        basis, pairinteraction::OperatorType::ELECTRIC_DIPOLE, 0);
    pairinteraction::OperatorAtom<std::complex<float>> dipole_p(
        basis, pairinteraction::OperatorType::ELECTRIC_DIPOLE, 1);
    pairinteraction::OperatorAtom<std::complex<float>> dipole_m(
        basis, pairinteraction::OperatorType::ELECTRIC_DIPOLE, -1);

    if (dipole_0.get_matrix().rows() != 64) {
        SPDLOG_ERROR("Wrong dimension.");
        return 1;
    }

    if (dipole_p.get_matrix().rows() != 64) {
        SPDLOG_ERROR("Wrong dimension.");
        return 1;
    }

    if (dipole_m.get_matrix().rows() != 64) {
        SPDLOG_ERROR("Wrong dimension.");
        return 1;
    }

    if (dipole_0.get_matrix().nonZeros() != 288) {
        SPDLOG_ERROR("Wrong number of non-zeros.");
        return 1;
    }

    if (dipole_p.get_matrix().nonZeros() != 288) {
        SPDLOG_ERROR("Wrong number of non-zeros.");
        return 1;
    }

    if (dipole_m.get_matrix().nonZeros() != 288) {
        SPDLOG_ERROR("Wrong number of non-zeros.");
        return 1;
    }

    return 0;
}
