par = {}

def push_element(key, value, name):
    par[key] = {
        'value' : value,
        'name'  : name,
        'type'  : type(value),
        'wx'    : None,
    }

push_element('species1', 'Rb', 'Element 1')
push_element('n1', 66, 'n1')
push_element('l1', 0, 'l1')
push_element('j1', 0.5, 'j1')
push_element('m1', 0.5, 'm1')
push_element('species2', 'Rb', 'Element 2')
push_element('n2', 66, 'n2')
push_element('l2', 0, 'l2')
push_element('j2', 0.5, 'j2')
push_element('m2', -0.5, 'm2')

push_element('steps', 21, 'Simulation steps')
push_element('minR', 0, 'Minimal distance (au)')
push_element('maxR', 1e5, 'Maximal distance (au)')
push_element('minBx', 0, 'Minimal b-field in X-direction (au)')
push_element('minBy', 0, 'Minimal b-field in Y-direction (au)')
push_element('minBz', 2.127191e-8, 'Minimal b-field in Z-direction (au)')
push_element('maxBx', 0, 'Maximal b-field in X-direction (au)')
push_element('maxBy', 0, 'Maximal b-field in Y-direction (au)')
push_element('maxBz', 2.127191e-8, 'Maximal b-field in Z-direction (au)')
push_element('minEx', 0, 'Minimal e-field in X-direction (au)')
push_element('minEy', 0, 'Minimal e-field in Y-direction (au)')
push_element('minEz', 0, 'Minimal e-field in Z-direction (au)')
push_element('maxEx', 0, 'Maximal e-field in X-direction (au)')
push_element('maxEy', 0, 'Maximal e-field in Y-direction (au)')
push_element('maxEz', 1e-11, 'Maximal e-field in Z-direction (au)')

push_element('deltaN', 4, 'ΔN')
push_element('deltaL', 10, 'ΔL')
push_element('deltaJ', 100, 'ΔJ')
push_element('deltaM', 1, 'ΔM')
push_element('deltaE', 0.7e-5, 'ΔEnergy (au)')

'''push_element('element', 'rb87', 'Element')
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
push_element('slp_enablesquaregrid', True, 'ScaLAPACK Square Grid')'''
