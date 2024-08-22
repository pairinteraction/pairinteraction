#include "./Database.py.hpp"

#include "pairinteraction/database/AtomDescriptionByParameters.hpp"
#include "pairinteraction/database/AtomDescriptionByRanges.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"
#include "pairinteraction/utils/Range.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>

namespace nb = nanobind;
using namespace pairinteraction;

template <typename T>
static void declare_range(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "Range" + type_name;
    nb::class_<Range<T>> pyclass(m, pylass_name.c_str());
    pyclass.def(nb::init<>())
        .def(nb::init<T, T>())
        .def("min", &Range<T>::min)
        .def("max", &Range<T>::max)
        .def("is_finite", &Range<T>::is_finite);
}

template <typename T>
static void declare_atom_description_by_parameters(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "AtomDescriptionByParameters" + type_name;
    nb::class_<AtomDescriptionByParameters<T>> pyclass(m, pylass_name.c_str());
    pyclass.def(nb::init<>())
        .def_rw("energy", &AtomDescriptionByParameters<T>::energy)
        .def_rw("quantum_number_f", &AtomDescriptionByParameters<T>::quantum_number_f)
        .def_rw("quantum_number_m", &AtomDescriptionByParameters<T>::quantum_number_m)
        .def_rw("parity", &AtomDescriptionByParameters<T>::parity)
        .def_rw("quantum_number_n", &AtomDescriptionByParameters<T>::quantum_number_n)
        .def_rw("quantum_number_nu", &AtomDescriptionByParameters<T>::quantum_number_nu)
        .def_rw("quantum_number_l", &AtomDescriptionByParameters<T>::quantum_number_l)
        .def_rw("quantum_number_s", &AtomDescriptionByParameters<T>::quantum_number_s)
        .def_rw("quantum_number_j", &AtomDescriptionByParameters<T>::quantum_number_j);
}

template <typename T>
static void declare_atom_description_by_ranges(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "AtomDescriptionByRanges" + type_name;
    nb::class_<AtomDescriptionByRanges<T>> pyclass(m, pylass_name.c_str());
    pyclass.def(nb::init<>())
        .def_rw("parity", &AtomDescriptionByRanges<T>::parity)
        .def_rw("range_energy", &AtomDescriptionByRanges<T>::range_energy)
        .def_rw("range_quantum_number_f", &AtomDescriptionByRanges<T>::range_quantum_number_f)
        .def_rw("range_quantum_number_m", &AtomDescriptionByRanges<T>::range_quantum_number_m)
        .def_rw("range_quantum_number_n", &AtomDescriptionByRanges<T>::range_quantum_number_n)
        .def_rw("range_quantum_number_nu", &AtomDescriptionByRanges<T>::range_quantum_number_nu)
        .def_rw("range_quantum_number_l", &AtomDescriptionByRanges<T>::range_quantum_number_l)
        .def_rw("range_quantum_number_s", &AtomDescriptionByRanges<T>::range_quantum_number_s)
        .def_rw("range_quantum_number_j", &AtomDescriptionByRanges<T>::range_quantum_number_j);
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
        .def(nb::init<bool>())
        .def(nb::init<std::filesystem::path>())
        .def(nb::init<bool, bool, std::filesystem::path>())
        .def("get_availability_of_species", &Database::get_availability_of_species)
        .def("get_availability_of_wigner_table", &Database::get_availability_of_wigner_table)
        .def_static("get_global_instance", nb::overload_cast<>(&Database::get_global_instance))
        .def_static("get_global_instance", nb::overload_cast<bool>(&Database::get_global_instance))
        .def_static("get_global_instance",
                    nb::overload_cast<std::filesystem::path>(&Database::get_global_instance))
        .def_static("get_global_instance",
                    nb::overload_cast<bool, bool, std::filesystem::path>(&Database::get_global_instance))
        .def("get_ket",
             nb::overload_cast<std::string, const AtomDescriptionByParameters<float> &>(
                 &Database::get_ket<float>))
        .def("get_ket",
             nb::overload_cast<std::string, const AtomDescriptionByParameters<double> &>(
                 &Database::get_ket<double>))
        .def("get_basis",
             nb::overload_cast<std::string, const AtomDescriptionByRanges<float> &,
                               std::vector<size_t>>(&Database::get_basis<float>))
        .def("get_basis",
             nb::overload_cast<std::string, const AtomDescriptionByRanges<double> &,
                               std::vector<size_t>>(&Database::get_basis<double>))
        .def("get_basis",
             nb::overload_cast<std::string, const AtomDescriptionByRanges<float> &,
                               std::vector<size_t>>(&Database::get_basis<std::complex<float>>))
        .def("get_basis",
             nb::overload_cast<std::string, const AtomDescriptionByRanges<double> &,
                               std::vector<size_t>>(&Database::get_basis<std::complex<double>>))
        .def("get_operator",
             nb::overload_cast<std::shared_ptr<const BasisAtom<float>>, OperatorType, int>(
                 &Database::get_operator<float>))
        .def("get_operator",
             nb::overload_cast<std::shared_ptr<const BasisAtom<double>>, OperatorType, int>(
                 &Database::get_operator<double>))
        .def("get_operator",
             nb::overload_cast<std::shared_ptr<const BasisAtom<std::complex<float>>>, OperatorType,
                               int>(&Database::get_operator<std::complex<float>>))
        .def("get_operator",
             nb::overload_cast<std::shared_ptr<const BasisAtom<std::complex<double>>>, OperatorType,
                               int>(&Database::get_operator<std::complex<double>>));
}

void bind_database(nb::module_ &m) {
    declare_range<int>(m, "Int");
    declare_range<float>(m, "Float");
    declare_range<double>(m, "Double");
    declare_atom_description_by_parameters<float>(m, "Float");
    declare_atom_description_by_parameters<double>(m, "Double");
    declare_atom_description_by_ranges<float>(m, "Float");
    declare_atom_description_by_ranges<double>(m, "Double");
    declare_availability_species(m);
    declare_availability_wigner(m);
    declare_database(m);
}
