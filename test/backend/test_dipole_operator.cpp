#include "pintr/pintr.hpp"
#include "pintr/utils/args.hpp"

#include <filesystem>
#include <spdlog/spdlog.h>

int main(int argc, char **argv) {
    // Call the setup function to configure logging
    pintr::setup();

    // Create a database instance
    std::filesystem::path databasedir;
    bool download_missing = false;

    for (int i = 1; i < argc; ++i) {
        bool found = pintr::args::parse_download_missing(i, argc, argv, download_missing);
        if (!found) {
            pintr::args::parse_database(i, argc, argv, databasedir);
        }
    }

    pintr::Database database(download_missing, true, databasedir);

    // Create a dipole operator coupling two specific states
    auto ket1 = pintr::KetAtomCreator<float>()
                    .set_species("Rb")
                    .set_quantum_number_n(60)
                    .set_quantum_number_l(0)
                    .set_quantum_number_j(0.5)
                    .set_quantum_number_m(0.5)
                    .create(database);

    auto ket2 = pintr::KetAtomCreator<float>()
                    .set_species("Rb")
                    .set_quantum_number_n(60)
                    .set_quantum_number_l(1)
                    .set_quantum_number_j(0.5)
                    .set_quantum_number_m(0.5)
                    .create(database);

    auto basis_ket1_ket2 =
        pintr::BasisAtomCreator<float>().add_ket(ket1).add_ket(ket2).create(database);

    pintr::OperatorAtom<float> dipole_ket1_ket2(basis_ket1_ket2,
                                                pintr::OperatorType::ELECTRIC_DIPOLE, 0);
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
    auto basis = pintr::BasisAtomCreator<std::complex<float>>()
                     .set_species("Sr88_singlet")
                     .restrict_quantum_number_n(60, 63)
                     .restrict_quantum_number_l(0, 3)
                     .create(database);

    pintr::OperatorAtom<std::complex<float>> dipole_0(basis, pintr::OperatorType::ELECTRIC_DIPOLE,
                                                      0);
    pintr::OperatorAtom<std::complex<float>> dipole_p(basis, pintr::OperatorType::ELECTRIC_DIPOLE,
                                                      1);
    pintr::OperatorAtom<std::complex<float>> dipole_m(basis, pintr::OperatorType::ELECTRIC_DIPOLE,
                                                      -1);

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
