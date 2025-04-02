// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

namespace pairinteraction {
enum class OperatorType {
    ZERO,                     // Empty operator
    ENERGY,                   // Unperturbed Hamiltonian
    ELECTRIC_DIPOLE,          // Dipole operator
    ELECTRIC_QUADRUPOLE,      // Quadrupole operator
    ELECTRIC_QUADRUPOLE_ZERO, // Part of the diamagnetic operator and quadrupole near surfaces
    ELECTRIC_OCTUPOLE,        // Octupole operator
    MAGNETIC_DIPOLE,          // Magnetic dipole operator
    IDENTITY,                 // Identity operator
    ARBITRARY                 // Arbitrary operator
};
} // namespace pairinteraction
