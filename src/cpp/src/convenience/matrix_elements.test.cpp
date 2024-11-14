#include "pairinteraction/convenience/matrix_elements.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerEigen.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/Range.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {

DOCTEST_TEST_CASE("Test the calculation of matrix elements") {
    auto &database = Database::get_global_instance();

    auto ket_s = KetAtomCreator<double>()
                     .set_species("Rb")
                     .set_quantum_number_n(60)
                     .set_quantum_number_l(0)
                     .set_quantum_number_j(0.5)
                     .set_quantum_number_m(0.5)
                     .create(database);

    auto ket_p = KetAtomCreator<double>()
                     .set_species("Rb")
                     .set_quantum_number_n(60)
                     .set_quantum_number_l(1)
                     .set_quantum_number_j(0.5)
                     .set_quantum_number_m(0.5)
                     .create(database);

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(59, 61)
                     .restrict_quantum_number_l(0, 1)
                     .restrict_quantum_number_m(0.5, 0.5)
                     .create(database);

    SystemAtom<double> system(basis);

    DOCTEST_SUBCASE("Test calculate_energy") {
        double energy1 = calculate_energy(ket_s);
        double energy2 = calculate_energy(ket_s, system);
        double energy3 = ket_s->get_energy();

        auto state = system.get_basis()->get_corresponding_state(ket_s);
        auto matrix_element =
            ket_s->get_database().get_matrix_elements(state, state, OperatorType::ENERGY, 0);
        DOCTEST_CHECK(matrix_element.rows() == 1);
        DOCTEST_CHECK(matrix_element.cols() == 1);
        double energy4 = matrix_element.coeff(0, 0);

        DOCTEST_CHECK(std::abs(energy1 - energy4) < 1e-11);
        DOCTEST_CHECK(std::abs(energy2 - energy4) < 1e-11);
        DOCTEST_CHECK(std::abs(energy3 - energy4) < 1e-11);
    }

    DOCTEST_SUBCASE("Test calculate_electric_dipole_matrix_element") {
        double dipole1 = calculate_electric_dipole_matrix_element(ket_s, ket_p, 0);
        double dipole2 = calculate_electric_dipole_matrix_element(ket_s, ket_p, system, 0);
        DOCTEST_CHECK(std::abs(dipole1 - dipole2) < 1e-11);
    }

    DOCTEST_SUBCASE("Test calculate_electric_dipole_matrix_element with an induced dipole") {
        double dipole = calculate_electric_dipole_matrix_element(ket_s, ket_s, system, 0);
        DOCTEST_CHECK(std::abs(dipole) < 1e-11);

        system.set_electric_field({0, 0, 1.9446903811524456e-10});
        system.diagonalize(DiagonalizerEigen<double>());
        dipole = calculate_electric_dipole_matrix_element(ket_s, ket_s, system, 0);
        DOCTEST_CHECK(std::abs(dipole) > 1e-11);
    }
}

} // namespace pairinteraction
