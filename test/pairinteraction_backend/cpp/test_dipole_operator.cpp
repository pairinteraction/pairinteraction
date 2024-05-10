#include "basis/BasisAtom.hpp"
#include "basis/BasisAtomCreator.hpp"
#include "database/Database.hpp"
#include "enums/OperatorType.hpp"
#include "ket/KetAtom.hpp"
#include "ket/KetAtomCreator.hpp"
#include "operator/OperatorAtom.hpp"
#include "setup.hpp"

#include <spdlog/spdlog.h>

int main() {
    // Call the setup function to configure logging
    setup();

    // Get a database instance
    Database &database = Database::get_global_instance();

    // Create a dipole operator coupling two specific states
    auto ket1 = KetAtomCreator<float>()
                    .set_species("Rb")
                    .set_quantum_number_n(60)
                    .set_quantum_number_l(0)
                    .set_quantum_number_j(0.5)
                    .set_quantum_number_m(0.5)
                    .create(database);

    auto ket2 = KetAtomCreator<float>()
                    .set_species("Rb")
                    .set_quantum_number_n(60)
                    .set_quantum_number_l(1)
                    .set_quantum_number_j(0.5)
                    .set_quantum_number_m(0.5)
                    .create(database);

    auto basis_ket1_ket2 = BasisAtomCreator<float>().add_ket(ket1).add_ket(ket2).create(database);

    OperatorAtom<float> dipole_ket1_ket2(basis_ket1_ket2, OperatorType::ELECTRIC_DIPOLE, 0);
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
    auto basis = BasisAtomCreator<std::complex<float>>()
                     .set_species("Sr88_singlet")
                     .restrict_quantum_number_n(60, 63)
                     .restrict_quantum_number_l(0, 3)
                     .create(database);

    OperatorAtom<std::complex<float>> dipole_0(basis, OperatorType::ELECTRIC_DIPOLE, 0);
    OperatorAtom<std::complex<float>> dipole_p(basis, OperatorType::ELECTRIC_DIPOLE, 1);
    OperatorAtom<std::complex<float>> dipole_m(basis, OperatorType::ELECTRIC_DIPOLE, -1);

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
