// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./LoggerBridge.hpp"
#include "./basis/Basis.py.hpp"
#include "./database/Database.py.hpp"
#include "./diagonalize/Diagonalizer.py.hpp"
#include "./enums/FloatType.py.hpp"
#include "./enums/OperatorType.py.hpp"
#include "./enums/Parity.py.hpp"
#include "./enums/TransformationType.py.hpp"
#include "./interfaces/DiagonalizerInterface.py.hpp"
#include "./interfaces/TransformationBuilderInterface.py.hpp"
#include "./ket/Ket.py.hpp"
#include "./operator/Operator.py.hpp"
#include "./paths.py.hpp"
#include "./system/System.py.hpp"
#include "./tools/run_unit_tests.py.hpp"
#include "./version.py.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

NB_MODULE(_backend, m) // NOLINT
{
    // https://nanobind.readthedocs.io/en/latest/faq.html#why-am-i-getting-errors-about-leaked-functions-and-types
    nb::set_leak_warnings(false);

    // wrap the get_pending_logs method of the logger bridge instance
    static LoggerBridge bridge;
    nb::class_<LoggerBridge::LogEntry>(m, "LogEntry")
        .def_ro("level", &LoggerBridge::LogEntry::level)
        .def_ro("message", &LoggerBridge::LogEntry::message);
    m.def("get_pending_logs", []() { return bridge.get_pending_logs(); });

    // enums
    bind_operator_type(m);
    bind_parity(m);
    bind_transformation_type(m);
    bind_float_type(m);

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

    // tools
    bind_run_unit_tests(m);

    // paths
    bind_paths(m);

    // version
    bind_version(m);
}
