#include "./basis/Basis.py.hpp"
#include "./database/Database.py.hpp"
#include "./enums/OperatorType.py.hpp"
#include "./enums/TransformationType.py.hpp"
#include "./interfaces/TransformationBuilderInterface.py.hpp"
#include "./ket/Ket.py.hpp"
#include "./operator/Operator.py.hpp"
#include "./system/System.py.hpp"
#include "pairinteraction/tools/setup.hpp"

#include <nanobind/nanobind.h>

NB_MODULE(backend, m) // NOLINT
{
    pairinteraction::setup();

    // enums
    bind_operator_type(m);
    bind_transformation_type(m);

    // interfaces
    bind_transformation_builder_interface(m);

    // operator
    bind_operator(m);

    // database
    bind_database(m);

    // ket
    bind_ket(m);

    // basis
    bind_basis(m);

    // system
    bind_system(m);
}
