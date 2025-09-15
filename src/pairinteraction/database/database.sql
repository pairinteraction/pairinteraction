-- Version 2.0

-- Unit System: Atomic Units
--
-- Base units:
-- - length: Bohr radius
-- - time: reduced Planck constant divided by Hartree energy
-- - temperature: Hartree energy divided by Boltzmann constant
-- - mass: electron mass
-- - current: elementary charge times Hartree energy divided by reduced Planck constant
--
-- Fundamental constants in this unit system:
-- - reduced Planck constant:   1
-- - Coulomb constant:          1
-- - Bohr magneton:             0.5
-- - elementary charge:         1
-- - electron mass:             1
--

BEGIN TRANSACTION;

-- Create a table of states
-- The energies are defined as E=E_i-Ry/nu_i^2, where E_i is the channel-dependent ionization threshold and nu_i is the
-- effective principal quantum number of the channel i (note that by using this definition, the energies of Rydberg
-- states are positive).
-- The effective principal quantum number of the state is defined as nu=sqrt(Ry/(I-E)). In the case of a single channel,
-- "I" is the ionization threshold of that channel and nu=n-delta where delta is the quantum defect. In the case of multiple
-- channels, "I" is the ionization threshold of the channel that corresponds to the core being in the ionic ground state and
-- has maximum total angular momentum f_core (this ensures that for stretched states with maximum F that consist only of a
-- single channel, the energy agrees with the energy obtainable by single-channel quantum defect theory).
-- The principal quantum number n is, in principle, not well defined for states that consist of multiple channels since
-- these states comprise superpositions of different radial wave functions. Instead, we state the principal quantum number
-- of the most dominant channel.
CREATE TABLE states (
    id INTEGER PRIMARY KEY,
    energy REAL NOT NULL, -- the energy of the state: E=I-Ry/nu^2
    parity INTEGER NOT NULL, -- the parity of the state
    n INTEGER NOT NULL, -- the principal quantum number n
    nu REAL NOT NULL, -- the effective principal quantum number nu
    f REAL NOT NULL, -- the total momentum quantum number f (int or half-int)
    exp_nui REAL NOT NULL, -- the expectation value of the effective principal quantum numbers nu_i of the channels
    exp_l REAL NOT NULL, -- the expectation value of the orbital quantum number l of all valence electrons
    exp_j REAL NOT NULL, -- the expectation value of the total angular quantum number j of all valence electrons
    exp_s REAL NOT NULL, -- the expectation value of the total spin quantum number s of all valence electrons
    exp_l_ryd REAL NOT NULL, -- the expectation value of the orbital quantum number l_{Ryd} of the Rydberg electron
    exp_j_ryd REAL NOT NULL, -- the expectation value of the total angular quantum number j_{Ryd} of the Rydberg electron
    std_nui REAL NOT NULL, -- the standard deviation of the effective principal quantum numbers nu_i of the channels
    std_l REAL NOT NULL, -- the standard deviation of the orbital quantum number l of all valence electrons
    std_j REAL NOT NULL, -- the standard deviation of the total angular quantum number j of all valence electrons
    std_s REAL NOT NULL, -- the standard deviation of the total spin quantum number s of all valence electrons
    std_l_ryd REAL NOT NULL, -- the standard deviation of the orbital quantum number l_{Ryd} of the Rydberg electron
    std_j_ryd REAL NOT NULL, -- the standard deviation of the total angular quantum number j_{Ryd} of the Rydberg electron
    is_j_total_momentum BOOLEAN NOT NULL, -- whether j is the total momentum quantum number, otherwise f is used for printing the state
    is_calculated_with_mqdt BOOLEAN NOT NULL, -- whether multi-channel quantum defect theory was used
    underspecified_channel_contribution REAL NOT NULL -- the contribution of channels whose quantum numbers are not exactly known
);

CREATE TRIGGER check_energy
BEFORE INSERT ON states
BEGIN
    SELECT RAISE(FAIL, "id and energy must be inserted in increasing order")
    FROM states
    WHERE id > NEW.id OR energy > NEW.energy;
END;

-- Create a table of wigner 3j symbols multiplied by the phase factor (-1)^(f-m)
CREATE TABLE wigner (
    f_initial REAL, -- the total momentum quantum number of the initial state
    f_final REAL, -- the total momentum quantum number of the final state
    m_initial REAL,
    m_final REAL,
    kappa INTEGER,
    q INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (f_initial, f_final, m_initial, m_final, kappa, q)
);

-- Create a table of matrix elements of the reduced electric dipole operator e*r*sqrt(4*pi/3)*Y_1
CREATE TABLE matrix_elements_d (
    id_initial INTEGER,
    id_final INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (id_initial, id_final)
);

-- Create a table of matrix elements of the reduced electric quadrupole operator e*r^2*sqrt(4*pi/5)*Y_2
CREATE TABLE matrix_elements_q (
    id_initial INTEGER,
    id_final INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (id_initial, id_final)
);

-- Create a table of matrix elements of the reduced operator e*r^2*sqrt(4*pi)*Y_0, occurring in
-- the diamagnetic interaction and electric quadrupole interaction near surfaces
CREATE TABLE matrix_elements_q0 (
    id_initial INTEGER,
    id_final INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (id_initial, id_final)
);

-- Create a table of matrix elements of the reduced electric octupole operator e*r^3*sqrt(4*pi/7)*Y_3
CREATE TABLE matrix_elements_o (
    id_initial INTEGER,
    id_final INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (id_initial, id_final)
);

-- Create a table of matrix elements of the reduced magnetic dipole operator
-- \mu = - \mu_B / \hbar * (g_L * L + g_S * S) + \mu_N / \hbar g_I I,
-- taking into account all spins: Rydberg electron, core electron, and nucleus
CREATE TABLE matrix_elements_mu (
    id_initial INTEGER,
    id_final INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (id_initial, id_final)
);

-- Create indices to speed up queries
CREATE INDEX idx_states_energy ON states (energy);
CREATE INDEX idx_wigner_f_initial ON wigner (f_initial);
CREATE INDEX idx_wigner_f_final ON wigner (f_final);
CREATE INDEX idx_matrix_elements_d_id_initial ON matrix_elements_d (id_initial);
CREATE INDEX idx_matrix_elements_d_id_final ON matrix_elements_d (id_final);
CREATE INDEX idx_matrix_elements_q_id_initial ON matrix_elements_q (id_initial);
CREATE INDEX idx_matrix_elements_q_id_final ON matrix_elements_q (id_final);
CREATE INDEX idx_matrix_elements_q0_id_initial ON matrix_elements_q0 (id_initial);
CREATE INDEX idx_matrix_elements_q0_id_final ON matrix_elements_q0 (id_final);
CREATE INDEX idx_matrix_elements_o_id_initial ON matrix_elements_o (id_initial);
CREATE INDEX idx_matrix_elements_o_id_final ON matrix_elements_o (id_final);
CREATE INDEX idx_matrix_elements_mu_id_initial ON matrix_elements_mu (id_initial);
CREATE INDEX idx_matrix_elements_mu_id_final ON matrix_elements_mu (id_final);

COMMIT;
