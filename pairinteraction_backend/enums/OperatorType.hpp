#pragma once

enum class OperatorType {
    ENERGY,         // Unperturbed Hamiltonian
    DIPOLE,         // Dipole operator
    QUADRUPOLE,     // Quadrupole operator
    OCTUPOLE,       // Octupole operator
    MAGNETICDIPOLE, // Magnetic dipole operator
    DIAMAGNETIC     // Part of the diamagnetic operator
};
