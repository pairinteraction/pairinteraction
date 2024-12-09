#include "./basis/Basis.py.hpp"
#include "./convenience/matrix_elements.py.hpp"
#include "./database/Database.py.hpp"
#include "./diagonalizer/Diagonalizer.py.hpp"
#include "./enums/OperatorType.py.hpp"
#include "./enums/Parity.py.hpp"
#include "./enums/TransformationType.py.hpp"
#include "./interfaces/DiagonalizerInterface.py.hpp"
#include "./interfaces/TransformationBuilderInterface.py.hpp"
#include "./ket/Ket.py.hpp"
#include "./operator/Operator.py.hpp"
#include "./system/System.py.hpp"
#include "./test.py.hpp"
#include "pairinteraction/tools/setup.hpp"

#include <nanobind/nanobind.h>

NB_MODULE(_backend, m) // NOLINT
{
    // https://nanobind.readthedocs.io/en/latest/faq.html#why-am-i-getting-errors-about-leaked-functions-and-types
    nanobind::set_leak_warnings(false);

    pairinteraction::setup();

    // enums
    bind_operator_type(m);
    bind_parity(m);
    bind_transformation_type(m);

    // interfaces
    bind_diagonalizer_interface(m);
    bind_transformation_builder_interface(m);

    // operator
    bind_operator(m);

    // database
    bind_database(m);

    // diagonalizer
    bind_diagonalizer(m);

    // ket
    bind_ket(m);

    // basis
    bind_basis(m);

    // system
    bind_system(m);

    // test
    bind_test(m);

    // convenience functions
    bind_matrix_elements(m);
}
