#!/usr/bin/env python3

import numpy as np
from scipy import sparse

# --- methods to read energies and basis from files ---

typeIds = {1008: 'int8', 1016 : 'int16', 1032 : 'int32', 1064 : 'int64', 1108 : 'uint8', 1116 : 'uint16', 1132 : 'uint32', \
    1164 : 'int64', 2032 : 'float32', 2064 : 'float64'}
type_t = 'uint16'
csr_not_csc = 0x01; # xxx0: csc, xxx1: csr
complex_not_real = 0x02; # xx0x: real, xx1x: complex

def readNumber(f, sz = 1):
    datatype = typeIds[np.fromfile(f, dtype=np.dtype(type_t), count=1)[0]]
    return np.squeeze(np.fromfile(f, dtype=np.dtype(datatype), count=sz))

def readVector(f):
    size = readNumber(f)
    return readNumber(f, size)

def readMatrix(f):
    flags = readNumber(f)
    rows = readNumber(f)
    cols = readNumber(f)
    if flags & complex_not_real: data = readVector(f) + readVector(f)*1j
    else: data = readVector(f)
    indices = readVector(f)
    indptr = np.append(readVector(f),len(data))
    if flags & csr_not_csc: return sparse.csr_matrix((data, indices, indptr), shape=(rows, cols))
    else: return sparse.csc_matrix((data, indices, indptr), shape=(rows, cols))

def load(filename):
    with open(filename,'rb') as f:
        energies = readMatrix(f).diagonal()
        basis = readMatrix(f)
        return energies, basis

def load2(filename):
    with open(filename,'rb') as f:
        energies = readMatrix(f)
        return energies


# --- usage example ---

# create output directory if necessary
import os
folder = "../build/output"
if not os.path.exists(folder):
    os.makedirs(folder)

"""# delete old files (otherwise the c++ program uses them, cache)
for the_file in os.listdir(folder):
    file_path = os.path.join(folder, the_file)
    if os.path.isfile(file_path): os.unlink(file_path)

# run the program
import subprocess
p = subprocess.Popen(["make","check"], stdout=subprocess.PIPE)
for line in iter(p.stdout.readline, b""):
    print (line.decode('utf-8'), end="")
    if not line: break"""
    

# plot lines
path_lines = os.path.join(folder,"lines.mat")
linessparse = load2(path_lines).T* 6579683.920729 #au to MHz
nSteps = 7*4
field = np.arange(nSteps)/(nSteps-1)*1e-11*5.14220652e11/100 #au to V/cm

lines = linessparse.todense()
lines_mask = np.ones_like(lines,dtype=np.bool)
lines_mask[linessparse.nonzero()] = False
lines[lines_mask] = np.nan

import matplotlib.pyplot as plt
plt.plot(field,lines,'-')
plt.xlim(-0.005,0.05)
#plt.ylim(-2,2)
plt.ylim(2.0,4.8)
plt.show()




exit()

# load results
nSteps = 7*4
files = os.path.join(folder,"hamiltonian_one_{}.mat")
energies_sym = [None]*nSteps
energies_asym = [None]*nSteps
for step in range(nSteps):
    energies_sym[step] = load(files.format(step))[0].real

energies_sym = np.array(energies_sym)
energies_sym *= 6579683.920729 #au to MHz

field = np.arange(nSteps)/(nSteps-1)*1e-11*5.14220652e11/100 #au to V/cm


# plot results
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
plt.plot(field,energies_sym,'r-')
plt.plot(field,energies_asym,'b-')
#plt.ylim(-0.055,0.021)
plt.xlim(-0.005,0.05)
#plt.xlim(0,0.4)
#plt.ylim(-1,5)
plt.ylim(2.0,4.8)
#plt.ylim(2.65,3.95)
#plt.xlim(-0.0005,0.047)
plt.show()
#plt.savefig('test.png')
