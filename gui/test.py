#!/usr/bin/env python3

import numpy as np
from scipy import sparse

# --- methods to read energies and basis from files ---

typeIds = {1016 : 'int16', 1032 : 'int32', 1064 : 'int64', 1116 : 'uint16', 1132 : 'uint32', \
    1164 : 'int64', 2032 : 'float32', 2064 : 'float64'}
type_t = 'uint16'

def readNumber(f, sz = 1):
    datatype = typeIds[np.fromfile(f, dtype=np.dtype(type_t), count=1)[0]]
    return np.squeeze(np.fromfile(f, dtype=np.dtype(datatype), count=sz))

def readVector(f):
    size = readNumber(f)
    return readNumber(f, size)

def readMatrix(f):
    mode = readNumber(f)
    rows = readNumber(f)
    cols = readNumber(f)
    data = readVector(f)
    indices = readVector(f)
    indptr = np.append(readVector(f),len(data))
    if mode == 0: return sparse.csc_matrix((data, indices, indptr), shape=(rows, cols))
    elif mode == 1: return sparse.csr_matrix((data, indices, indptr), shape=(rows, cols))

def load(filename):
    with open(filename,'rb') as f:
        energies = readMatrix(f).diagonal()
        basis = readMatrix(f)
        return energies, basis
    return None, None


# --- usage example ---

# create output directory if necessary
import os
folder = "output"
if not os.path.exists(folder):
    os.makedirs(folder)

# delete old files (otherwise the c++ program uses them, cache)
for the_file in os.listdir(folder):
    file_path = os.path.join(folder, the_file)
    if os.path.isfile(file_path): os.unlink(file_path)

# run the program
import subprocess
p = subprocess.Popen(["make","check"], stdout=subprocess.PIPE)
for line in iter(p.stdout.readline, b""):
    print (line.decode('utf-8'), end="")
    if not line: break

# load results
nSteps = 100
files = os.path.join(folder,"hamiltonian_two_{}_{}.mat")
energies_sym = [None]*nSteps
energies_asym = [None]*nSteps
for step in range(nSteps):
    energies_sym[step] = load(files.format(step,0))[0]
    energies_asym[step] = load(files.format(step,1))[0]

# plot results
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.plot(energies_sym,'r-')
plt.plot(energies_asym,'b-')
plt.savefig('test.png')