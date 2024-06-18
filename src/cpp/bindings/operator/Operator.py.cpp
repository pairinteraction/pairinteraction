#include "Operator.py.hpp"

#include "Operator.hpp"
#include "OperatorAtom.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>

namespace nb = nanobind;

template <typename T>
static void declare_operator(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "Operator" + type_name;
    using basis_t = typename Operator<T>::basis_t;
    using scalar_t = typename Operator<T>::scalar_t;
    nb::class_<Operator<T>, TransformationBuilderInterface<scalar_t>> pyclass(m,
                                                                              pylass_name.c_str());
    pyclass.def(nb::init<std::shared_ptr<const basis_t>>())
        .def("get_basis", &Operator<T>::get_basis)
        .def("get_matrix", &Operator<T>::get_matrix)
        .def("get_transformation", &Operator<T>::get_transformation)
        .def("get_rotator", &Operator<T>::get_rotator)
        .def("get_sorter", &Operator<T>::get_sorter)
        .def("get_blocks", &Operator<T>::get_blocks)
        .def("transform",
             nb::overload_cast<const Transformation<scalar_t> &>(&Operator<T>::transformed,
                                                                 nb::const_))
        .def("transform", nb::overload_cast<const Sorting &>(&Operator<T>::transformed, nb::const_))
        .def(scalar_t() * nb::self)
        .def(nb::self * scalar_t())
        .def(nb::self / scalar_t())
        .def(nb::self + nb::self)
        .def(nb::self - nb::self);
}

template <typename T>
static void declare_operator_atom(nb::module_ &m, std::string type_name) {
    std::string pylass_name = "OperatorAtom" + type_name;
    using basis_t = typename OperatorAtom<T>::basis_t;
    nb::class_<OperatorAtom<T>, Operator<OperatorAtom<T>>> pyclass(m, pylass_name.c_str());
    pyclass.def(nb::init<std::shared_ptr<const basis_t>>())
        .def(nb::init<std::shared_ptr<const basis_t>, OperatorType, int>());
}

void bind_operator(nb::module_ &m) {
    declare_operator<OperatorAtom<float>>(m, "OperatorAtomFloat");
    declare_operator<OperatorAtom<double>>(m, "OperatorAtomDouble");
    declare_operator<OperatorAtom<std::complex<float>>>(m, "OperatorAtomComplexFloat");
    declare_operator<OperatorAtom<std::complex<double>>>(m, "OperatorAtomComplexDouble");
    declare_operator_atom<float>(m, "Float");
    declare_operator_atom<double>(m, "Double");
    declare_operator_atom<std::complex<float>>(m, "ComplexFloat");
    declare_operator_atom<std::complex<double>>(m, "ComplexDouble");
}
