# SPDX-FileCopyrightText: 2024 Sebastian Weber, Henri Menke, and contributors. All rights reserved.
# SPDX-License-Identifier: GPL-3.0-or-later
import json

import numpy as np
from scipy import sparse


class BinaryLoader:
    def __init__(self):
        # types
        self.typeIds = {
            1008: "int8",
            1016: "int16",
            1032: "int32",
            1064: "int64",
            1108: "uint8",
            1116: "uint16",
            1132: "uint32",
            1164: "int64",
            2032: "float32",
            2064: "float64",
        }
        self.type_t = "uint16"

        # bit masks
        self.csr_not_csc = 0x01  # xxx0: csc, xxx1: csr
        self.complex_not_real = 0x02  # xx0x: real, xx1x: complex

    def readNumber(self, f, sz=None):
        datatype = self.typeIds[np.fromfile(f, dtype=np.dtype(self.type_t), count=1)[0]]
        if sz is None:
            return np.fromfile(f, dtype=np.dtype(datatype), count=1)[0]
        else:
            return np.fromfile(f, dtype=np.dtype(datatype), count=sz)

    def readVector(self, f):
        size = self.readNumber(f)
        return self.readNumber(f, size)

    def readMatrix(self, f):
        flags = self.readNumber(f)
        rows = self.readNumber(f)
        cols = self.readNumber(f)
        if flags & self.complex_not_real:
            data = self.readVector(f) + self.readVector(f) * 1j
        else:
            data = self.readVector(f)
        indices = self.readVector(f)
        indptr = np.append(self.readVector(f), len(data))
        if flags & self.csr_not_csc:
            return sparse.csr_matrix((data, indices, indptr), shape=(rows, cols))
        else:
            return sparse.csc_matrix((data, indices, indptr), shape=(rows, cols))


class Eigensystem(BinaryLoader):
    def __init__(self, filename):
        super().__init__()

        self._filename = filename
        self._shift = 0

        self._params = None
        self._energies = None
        self._basis = None

    @property
    def params(self):
        if self._params is None:
            with open(self._filename + ".json") as f:
                self._params = json.load(f)
        return self._params

    @property
    def energies(self):
        if self._energies is None:
            with open(self._filename + ".mat", "rb") as f:
                tmp = self.readMatrix(f)
                if tmp.size:
                    self._energies = np.real(tmp.diagonal())
                else:
                    self._energies = np.array([])
                self._shift = f.tell()
        return self._energies

    @property
    def basis(self):
        if self._basis is None:
            with open(self._filename + ".mat", "rb") as f:
                if self._shift > 0:
                    f.seek(self._shift, 0)
                else:
                    tmp = self.readMatrix(f)
                    if tmp.size:
                        self._energies = np.real(tmp.diagonal())
                    else:
                        self._energies = np.array([])
                self._basis = self.readMatrix(f)
        return self._basis
