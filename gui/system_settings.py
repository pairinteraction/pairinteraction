par = {}

def push_element(key, value, name):
    par[key] = {
        'value' : value,
        'name'  : name,
        'type'  : type(value),
        'wx'    : None,
    }

push_element('element', 'rb87', 'Element')
push_element('n1', 80, 'n1')
push_element('n2', 80, 'n2')
push_element('l1', 2, 'l1')
push_element('l2', 2, 'l2')
push_element('j1', 2.5, 'j1')
push_element('j2', 2.5, 'j2')
push_element('m1', 2.5, 'm1')
push_element('m2', 2.5, 'm2')
push_element('symmetry', 's', 'Symmetry')
push_element('angle', 0, 'Angle')
push_element('enable_dipdip', True, 'Dipole-Dipole')
push_element('enable_dipquad', True, 'Dipole-Quadrupole')
push_element('enable_quadquad', True, 'Quadrupole-Quadrupole')
push_element('efield_strength', 0, 'Electric Field [V/cm]')
push_element('efield_increasing', False, 'Increase Electric Field')
push_element('bfield_strength', 0, 'Magnetic Field [G]')
push_element('bfield_increasing', False, 'Increase Magnetic Field')
push_element('delta_n', 2, 'Δn')
push_element('delta_l', 2, 'Δl')
push_element('delta_m', 100, 'Δm')
push_element('delta_energy', 8, 'Δ energy [GHz]')
push_element('preserve_M', True, 'Preserve M')
push_element('preserve_submatrix', True, 'Preserve Submatrix')
push_element('preserve_parityL', False, 'Preserve Parity L')
#push_element('reducing_method', 'no', 'Reducing Method') # [not implemented yet]
#push_element('reducing_param1', 0, 'Reducing Parameter 1') # [not implemented yet]
#push_element('reducing_param1', 0, 'Reducing Parameter 1') # [not implemented yet]
push_element('distance_min', 2, 'Minimal Distance [µm]')
push_element('distance_max', 22, 'Maximal Distance [µm]')
push_element('distance_steps', 200, 'Distance Steps')
push_element('energyrange', 0.2, 'Energy Range [GHz]')
push_element('diagonalizer', 'ScaLAPACK', 'Diagonalizer')
push_element('cutoff_eigenVec', 1e-2, 'Cut-off Eigenvectors')
push_element('cutoff_interactionMat', 1e-12, 'Cut-off Interaction Matrix')
push_element('slepc_numeigenpairs', 100, 'SLEPc Eigenpairs')
push_element('slepc_solvertollerance', 1e-2, 'SLEPc Tolerance')
push_element('slepc_maxiteration', 500, 'SLEPc Maximal Iterations')
push_element('slp_blockdimension', 64, 'ScaLAPACK Block Dimension')
push_element('slp_enablesquaregrid', True, 'ScaLAPACK Square Grid')
