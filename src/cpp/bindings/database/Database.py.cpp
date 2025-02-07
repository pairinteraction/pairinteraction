#include "./Database.py.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/AtomDescriptionByParameters.hpp"
#include "pairinteraction/database/AtomDescriptionByRanges.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/utils/Range.hpp"

#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace pairinteraction;

template <typename T>
static void declare_range(nb::module_ &m, std::string const &type_name) {
    std::string pyclass_name = "Range" + type_name;
    nb::class_<Range<T>> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def(nb::init<T, T>())
        .def("min", &Range<T>::min)
        .def("max", &Range<T>::max)
        .def("is_finite", &Range<T>::is_finite);
}

static void declare_atom_description_by_parameters(nb::module_ &m) {
    std::string pyclass_name = "AtomDescriptionByParameters";
    nb::class_<AtomDescriptionByParameters> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def_rw("energy", &AtomDescriptionByParameters::energy)
        .def_rw("quantum_number_f", &AtomDescriptionByParameters::quantum_number_f)
        .def_rw("quantum_number_m", &AtomDescriptionByParameters::quantum_number_m)
        .def_rw("parity", &AtomDescriptionByParameters::parity)
        .def_rw("quantum_number_n", &AtomDescriptionByParameters::quantum_number_n)
        .def_rw("quantum_number_nu", &AtomDescriptionByParameters::quantum_number_nu)
        .def_rw("quantum_number_l", &AtomDescriptionByParameters::quantum_number_l)
        .def_rw("quantum_number_s", &AtomDescriptionByParameters::quantum_number_s)
        .def_rw("quantum_number_j", &AtomDescriptionByParameters::quantum_number_j);
}

static void declare_atom_description_by_ranges(nb::module_ &m) {
    std::string pyclass_name = "AtomDescriptionByRanges";
    nb::class_<AtomDescriptionByRanges> pyclass(m, pyclass_name.c_str());
    pyclass.def(nb::init<>())
        .def_rw("parity", &AtomDescriptionByRanges::parity)
        .def_rw("range_energy", &AtomDescriptionByRanges::range_energy)
        .def_rw("range_quantum_number_f", &AtomDescriptionByRanges::range_quantum_number_f)
        .def_rw("range_quantum_number_m", &AtomDescriptionByRanges::range_quantum_number_m)
        .def_rw("range_quantum_number_n", &AtomDescriptionByRanges::range_quantum_number_n)
        .def_rw("range_quantum_number_nu", &AtomDescriptionByRanges::range_quantum_number_nu)
        .def_rw("range_quantum_number_l", &AtomDescriptionByRanges::range_quantum_number_l)
        .def_rw("range_quantum_number_s", &AtomDescriptionByRanges::range_quantum_number_s)
        .def_rw("range_quantum_number_j", &AtomDescriptionByRanges::range_quantum_number_j);
}

static void declare_availability_species(nb::module_ &m) {
    nb::class_<Database::AvailabilitySpecies>(m, "DatabaseAvailabilitySpecies")
        .def(nb::init<std::string const &, bool, bool, bool>())
        .def_rw("name", &Database::AvailabilitySpecies::name)
        .def_rw("locally_available", &Database::AvailabilitySpecies::locally_available)
        .def_rw("up_to_date", &Database::AvailabilitySpecies::up_to_date)
        .def_rw("fully_downloaded", &Database::AvailabilitySpecies::fully_downloaded);
}

static void declare_availability_wigner(nb::module_ &m) {
    nb::class_<Database::AvailabilityWigner>(m, "DatabaseAvailabilityWigner")
        .def(nb::init<bool, bool>())
        .def_rw("locally_available", &Database::AvailabilityWigner::locally_available)
        .def_rw("up_to_date", &Database::AvailabilityWigner::up_to_date);
}

static void declare_database(nb::module_ &m) {
    nb::class_<Database>(m, "Database")
        .def(nb::init<>())
        .def(nb::init<bool>(), "download_missing"_a)
        .def(nb::init<std::filesystem::path>(), "database_dir"_a)
        .def(nb::init<bool, bool, std::filesystem::path>(), "download_missing"_a,
             "wigner_in_memory"_a, "database_dir"_a)
        .def("get_availability_of_species", &Database::get_availability_of_species)
        .def("get_availability_of_wigner_table", &Database::get_availability_of_wigner_table)
        .def(
            "get_ket",
            nb::overload_cast<std::string, const AtomDescriptionByParameters &>(&Database::get_ket))
        .def("get_basis",
             nb::overload_cast<std::string, const AtomDescriptionByRanges &, std::vector<size_t>>(
                 &Database::get_basis<double>))
        .def("get_basis",
             nb::overload_cast<std::string, const AtomDescriptionByRanges &, std::vector<size_t>>(
                 &Database::get_basis<std::complex<double>>))
        .def("get_matrix_elements",
             nb::overload_cast<std::shared_ptr<const BasisAtom<double>>,
                               std::shared_ptr<const BasisAtom<double>>, OperatorType, int>(
                 &Database::get_matrix_elements<double>))
        .def("get_matrix_elements",
             nb::overload_cast<std::shared_ptr<const BasisAtom<std::complex<double>>>,
                               std::shared_ptr<const BasisAtom<std::complex<double>>>, OperatorType,
                               int>(&Database::get_matrix_elements<std::complex<double>>));
}

void bind_database(nb::module_ &m) {
    declare_range<int>(m, "Int");
    declare_range<double>(m, "Double");
    declare_atom_description_by_parameters(m);
    declare_atom_description_by_ranges(m);
    declare_availability_species(m);
    declare_availability_wigner(m);
    declare_database(m);
}
