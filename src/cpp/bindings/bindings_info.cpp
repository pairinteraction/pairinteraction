// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./Info.py.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;

NB_MODULE(_info, m) // NOLINT
{
    nb::set_leak_warnings(false);
    bind_info(m);
}
