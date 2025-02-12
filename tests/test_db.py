"""Test receiving matrix elements from the databases."""

from pathlib import Path

import duckdb
import numpy as np
import pytest
from packaging.version import Version
from sympy.physics.wigner import wigner_3j

import pairinteraction.backend.real as pi
from tests.constants import GAUSS_IN_ATOMIC_UNITS, HARTREE_IN_GHZ, SPECIES_TO_NUCLEAR_SPIN, SUPPORTED_SPECIES

expected_operator = np.array([[3.58588117, 1.66420213], [1.66420213, 4.16645123]])


def fetch_id(n: int, l: float, f: float, s: float, connection: duckdb.duckdb.DuckDBPyConnection, table: str) -> int:
    return connection.execute(
        f"SELECT id FROM '{table}' WHERE n = {n} AND f = {f} ORDER BY (exp_s - {s})^2+(exp_l - {l})^2 LIMIT 1"
    ).fetchone()[0]


def fetch_wigner_element(
    f_initial: float,
    f_final: float,
    m_initial: float,
    m_final: float,
    kappa: int,
    q: int,
    connection: duckdb.duckdb.DuckDBPyConnection,
    table: str,
) -> float:
    result = connection.execute(
        f"SELECT val FROM '{table}' WHERE f_initial = {f_initial} AND f_final = {f_final} AND "
        f"m_initial = {m_initial} AND m_final = {m_final} AND kappa = {kappa} AND q = {q}"
    ).fetchone()
    return result[0] if result else 0


def fetch_reduced_matrix_element(
    id_initial: int, id_final: int, connection: duckdb.duckdb.DuckDBPyConnection, table: str
) -> float:
    result = connection.execute(
        f"SELECT val FROM '{table}' WHERE id_initial = {id_initial} AND id_final = {id_final}"
    ).fetchone()
    return result[0] if result else 0


@pytest.fixture(scope="module")
def connection() -> duckdb.duckdb.DuckDBPyConnection:
    with duckdb.connect(":memory:") as connection:
        yield connection


@pytest.mark.parametrize("swap_states", [False, True])
def test_database(connection: duckdb.duckdb.DuckDBPyConnection, swap_states: bool) -> None:
    """Test receiving matrix elements from the databases."""
    database = pi.Database.get_global_database()
    bfield_in_gauss = 1500

    # Define initial and final quantum states
    n_initial, n_final = 54, 54
    l_initial, l_final = 1, 1
    f_initial, f_final = 1, 0
    m_initial, m_final = 0, 0
    s_initial, s_final = 0.6, 1.0

    # Swap states if required by the test parameter
    if swap_states:
        n_initial, n_final = n_final, n_initial
        l_initial, l_final = l_final, l_initial
        f_initial, f_final = f_final, f_initial
        m_initial, m_final = m_final, m_initial
        s_initial, s_final = s_final, s_initial

    # Get the Zeeman interaction operator from the database using pairinteraction
    ket_initial = pi.KetAtom("Yb174_mqdt", n=n_initial, l=l_initial, f=f_initial, m=m_initial, s=s_initial)
    ket_final = pi.KetAtom("Yb174_mqdt", n=n_final, l=l_final, f=f_final, m=m_final, s=s_final)
    basis = pi.BasisAtom("Yb174_mqdt", additional_kets=[ket_initial, ket_final])
    operator = (
        pi.SystemAtom(basis)
        .set_magnetic_field([0, 0, bfield_in_gauss], unit="G")
        .enable_diamagnetism(True)
        .get_hamiltonian(unit="GHz")
    )
    operator -= np.diag(np.sort([ket_initial.get_energy(unit="GHz"), ket_final.get_energy(unit="GHz")]))
    assert np.allclose(operator, expected_operator)

    # Get the latest parquet files from the database directory
    parquet_files = {}
    parquet_versions = {}
    for path in list(Path(database.database_dir).glob("*.parquet")):
        name, version_str = path.stem.rsplit("_v", 1)
        version = Version(version_str)
        if name not in parquet_files or version > parquet_versions[name]:
            parquet_files[name] = path
            parquet_versions[name] = version
    assert "wigner" in parquet_files.keys()
    assert "Yb174_mqdt_states" in parquet_files.keys()
    assert "Yb174_mqdt_matrix_elements_mu" in parquet_files.keys()
    assert "Yb174_mqdt_matrix_elements_q" in parquet_files.keys()
    assert "Yb174_mqdt_matrix_elements_q0" in parquet_files.keys()

    # Obtain the ids of the initial and final states
    id_initial = fetch_id(n_initial, l_initial, f_initial, s_initial, connection, parquet_files["Yb174_mqdt_states"])
    assert id_initial == 362 if swap_states else 363

    id_final = fetch_id(n_final, l_final, f_final, s_final, connection, parquet_files["Yb174_mqdt_states"])
    assert id_final == 363 if swap_states else 362

    # Obtain a matrix element of the magnetic dipole operator (for the chosen kets, it is non-zero iff initial != final)
    kappa, q = 1, 0
    wigner_element = fetch_wigner_element(
        f_initial, f_final, m_initial, m_final, kappa, q, connection, parquet_files["wigner"]
    )
    assert np.isclose(
        wigner_element,
        float((-1) ** (f_final - m_final) * wigner_3j(f_final, kappa, f_initial, -m_final, q, m_initial)),
    )

    mu = fetch_reduced_matrix_element(id_initial, id_final, connection, parquet_files["Yb174_mqdt_matrix_elements_mu"])
    matrix_element = -wigner_element * mu * bfield_in_gauss * GAUSS_IN_ATOMIC_UNITS * HARTREE_IN_GHZ
    assert np.isclose(matrix_element, operator[0, 1])

    # Obtain a matrix element of the diamagnetic operator (for the chosen kets, it is non-zero iff initial == final)
    n_final = n_initial
    l_final = l_initial
    f_final = f_initial
    m_final = m_initial
    s_final = s_initial
    id_final = id_initial

    kappa, q = 0, 0
    wigner_element = fetch_wigner_element(
        f_initial, f_final, m_initial, m_final, kappa, q, connection, parquet_files["wigner"]
    )
    assert np.isclose(
        wigner_element,
        float((-1) ** (f_final - m_final) * wigner_3j(f_final, kappa, f_initial, -m_final, q, m_initial)),
    )

    q0 = fetch_reduced_matrix_element(id_initial, id_final, connection, parquet_files["Yb174_mqdt_matrix_elements_q0"])
    matrix_element = 1 / 12 * wigner_element * q0 * (bfield_in_gauss * GAUSS_IN_ATOMIC_UNITS) ** 2 * HARTREE_IN_GHZ

    kappa, q = 2, 0
    wigner_element = fetch_wigner_element(
        f_initial, f_final, m_initial, m_final, kappa, q, connection, parquet_files["wigner"]
    )
    assert np.isclose(
        wigner_element,
        float((-1) ** (f_final - m_final) * wigner_3j(f_final, kappa, f_initial, -m_final, q, m_initial)),
    )

    q = fetch_reduced_matrix_element(id_initial, id_final, connection, parquet_files["Yb174_mqdt_matrix_elements_q"])
    matrix_element -= 1 / 12 * wigner_element * q * (bfield_in_gauss * GAUSS_IN_ATOMIC_UNITS) ** 2 * HARTREE_IN_GHZ
    assert np.isclose(matrix_element, operator[0, 0] if swap_states else operator[1, 1])


@pytest.mark.parametrize("species", SUPPORTED_SPECIES)
def test_obtaining_kets(species: str) -> None:
    """Test obtaining kets from the database."""
    is_mqdt = species.endswith("_mqdt")
    is_single_valence_electron = species in ["Rb"]
    is_triplet = species in ["Sr88_triplet"]

    quantum_number_i = SPECIES_TO_NUCLEAR_SPIN[species] if is_mqdt else 0
    quantum_number_s = 0.5 if is_single_valence_electron else (1 if is_triplet else 0)
    quantum_number_f = quantum_number_i + quantum_number_s
    quantum_number_m = quantum_number_i + quantum_number_s

    # Obtain a ket from the database
    ket = pi.KetAtom(species, n=60, l=0, f=quantum_number_f, m=quantum_number_m, s=quantum_number_s)

    print(f"|{ket}>")

    # Check the result
    assert ket.species == species
    assert ket.n == 60 if not is_mqdt else abs(ket.n - 60) < 1
    assert ket.l == 0 if not is_mqdt else abs(ket.l - 0) < 1
    assert ket.f == quantum_number_f
    assert ket.m == quantum_number_m
    assert ket.s == quantum_number_s if not is_mqdt else abs(ket.s - quantum_number_s) < 1
    assert {
        "Rb": "Rb:60,S_1/2,1/2",
        "Sr88_singlet": "Sr88_singlet:60,S_0,0",
        "Sr88_triplet": "Sr88_triplet:60,S_1,1",
        "Sr88_mqdt": "Sr88:S=0.0,nu=56.7,L=0.0,J=0,0",
        "Sr87_mqdt": "Sr87:S=0.5,nu=56.6,L=0.0,F=9/2,9/2",
        "Yb174_mqdt": "Yb174:S=0.0,nu=55.7,L=0.0,J=0,0",
        "Yb171_mqdt": "Yb171:S=0.3,nu=55.7,L=0.0,F=1/2,1/2",
        "Yb173_mqdt": "Yb173:S=0.4,nu=56.4,L=0.0,F=5/2,5/2",
    }[species] == str(ket)
