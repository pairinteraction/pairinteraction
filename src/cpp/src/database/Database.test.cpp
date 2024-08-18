#include "pairinteraction/database/Database.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/AtomDescriptionByParameters.hpp"
#include "pairinteraction/database/AtomDescriptionByRanges.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/utils/streamed.hpp"

#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("get a KetAtom") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByParameters<float> description;
    description.quantum_number_n = 60;
    description.quantum_number_l = 0;
    description.quantum_number_m = 0;

    auto ket = database.get_ket<float>("Rb", description);

    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "KetAtom: {}", fmt::streamed(*ket));
}

DOCTEST_TEST_CASE("get a BasisAtom") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByRanges<float> description;
    description.range_quantum_number_n = {60, 60};
    description.range_quantum_number_l = {0, 1};

    auto basis = database.get_basis<float>("Rb", description, {});

    for (auto ket : *basis) {
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "KetAtom: {}", fmt::streamed(*ket));
    }
}

DOCTEST_TEST_CASE("get an OperatorAtom") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByRanges<float> description;
    description.range_quantum_number_n = {60, 60};
    description.range_quantum_number_l = {0, 1};

    auto basis = database.get_basis<float>("Rb", description, {});

    auto dipole = database.get_operator<float>(basis, OperatorType::ELECTRIC_DIPOLE, 0);

    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Number of basis states: {}",
                       basis->get_number_of_states());
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Number of non-zero entries: {}",
                       dipole.get_matrix().nonZeros());
}
} // namespace pairinteraction
