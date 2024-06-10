#include "TransformationBuilderInterface.py.hpp"

#include "TransformationBuilderInterface.hpp"

#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

template <typename T>
static void declare_transformation(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "Transformation" + type_name;
    nb::class_<Transformation<T>> pyclass(m, pylass_name.c_str());
    pyclass.def(nb::init<>())
        .def_rw("matrix", &Transformation<T>::matrix)
        .def_rw("transformation_type", &Transformation<T>::transformation_type);
}

static void declare_sorting(nb::module_ &m) {
    nb::class_<Sorting> pyclass(m, "Sorting");
    pyclass.def(nb::init<>())
        .def_rw("matrix", &Sorting::matrix)
        .def_rw("transformation_type", &Sorting::transformation_type);
}

static void declare_blocks(nb::module_ &m) {
    nb::class_<Blocks> pyclass(m, "Blocks");
    pyclass.def(nb::init<>())
        .def_rw("start", &Blocks::start)
        .def_rw("transformation_type", &Blocks::transformation_type);
}

template <typename T>
static void declare_transformation_builder_interface(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "TransformationBuilderInterface" + type_name;
    using real_t = typename TransformationBuilderInterface<T>::real_t;
    nb::class_<TransformationBuilderInterface<T>> pyclass(m, pylass_name.c_str());
    pyclass.def("get_rotator",
                nb::overload_cast<std::array<real_t, 3>, std::array<real_t, 3>>(
                    &TransformationBuilderInterface<T>::get_rotator, nb::const_));
}

void bind_transformation_builder_interface(nb::module_ &m) {
    declare_transformation<float>(m, "Float");
    declare_transformation<double>(m, "Double");
    declare_sorting(m);
    declare_blocks(m);
    declare_transformation_builder_interface<float>(m, "Float");
    declare_transformation_builder_interface<double>(m, "Double");
    declare_transformation_builder_interface<std::complex<float>>(m, "ComplexFloat");
    declare_transformation_builder_interface<std::complex<double>>(m, "ComplexDouble");
}
