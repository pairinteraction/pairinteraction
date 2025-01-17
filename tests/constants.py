"""
Constants used in the tests.
"""

SUPPORTED_SPECIES = [
    "Rb",
    "Sr88_singlet",
    "Sr88_triplet",
    "Sr88_mqdt",
    "Sr87_mqdt",
    "Yb174_mqdt",
    "Yb171_mqdt",
    "Yb173_mqdt",
]
SPECIES_TO_NUCLEAR_SPIN = {
    "Yb171_mqdt": 1 / 2,
    "Yb173_mqdt": 5 / 2,
    "Yb174_mqdt": 0,
    "Sr87_mqdt": 9 / 2,
    "Sr88_mqdt": 0,
}

HARTREE_IN_JOULES = 4.3597447222060e-18
HARTREE_IN_GHZ = 6579683.920502
HARTREE_IN_INVERSE_CM = 219474.63136320
VOLT_PER_CM_IN_ATOMIC_UNITS = 1 / 5.14220675112e9
GAUSS_IN_ATOMIC_UNITS = 1 / 2.35051757077e9
