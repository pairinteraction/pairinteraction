# SPDX-FileCopyrightText: 2024 Sebastian Weber, Henri Menke, and contributors. All rights reserved.
# SPDX-License-Identifier: GPL-3.0-or-later
import os
import pickle

import numpy as np

# === Calculate and cache Wigner d-matrix elements ===


class Wignerd:
    def __init__(self, cachedir) -> None:
        self.cachedir = cachedir
        self.wignerdict = {}
        self.cachupdatedict = {}

    def __del__(self) -> None:
        self.save()

    def save(self) -> None:
        for k, v in self.wignerdict.items():
            if not self.cachupdatedict[k]:
                continue
            path = os.path.join(self.cachedir, k)
            with open(path, "wb") as f:
                pickle.dump(v, f, pickle.HIGHEST_PROTOCOL)
                self.cachupdatedict[k] = False

    def calc(self, j, m2, m1, beta):
        sgn = 1

        if m1 + m2 < 0:
            m1 *= -1
            m2 *= -1
            sgn *= (-1) ** (m2 - m1)

        if m1 < m2:
            m1, m2 = m2, m1
            sgn *= (-1) ** (m2 - m1)

        bstring = f"{int(np.round(beta * 1e9)):+013d}.pkl"
        if bstring not in self.wignerdict.keys():
            path = os.path.join(self.cachedir, bstring)
            if os.path.exists(path):
                with open(path, "rb") as f:
                    self.wignerdict[bstring] = pickle.load(f)
            else:
                self.wignerdict[bstring] = {}
            self.cachupdatedict[bstring] = False

        mstring = f"{j}_{m1}_{m2}"
        if mstring not in self.wignerdict[bstring].keys():
            self.wignerdict[bstring][mstring] = float(np.real(self._eval_wignerd(j, m2, m1, beta)))
            self.cachupdatedict[bstring] = True

        return sgn * self.wignerdict[bstring][mstring]

    # This is a stripped down implementation of the WignerD symbols as
    # found in sympy.physics.quantum.spin
    #
    # Copyright (c) 2006-2016 SymPy Development Team
    #
    # All rights reserved.
    #
    # Redistribution and use in source and binary forms, with or without
    # modification, are permitted provided that the following conditions are met:
    #
    #   a. Redistributions of source code must retain the above copyright notice,
    #      this list of conditions and the following disclaimer.
    #   b. Redistributions in binary form must reproduce the above copyright
    #      notice, this list of conditions and the following disclaimer in the
    #      documentation and/or other materials provided with the distribution.
    #   c. Neither the name of SymPy nor the names of its contributors
    #      may be used to endorse or promote products derived from this software
    #      without specific prior written permission.
    #
    #
    # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    # AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    # ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
    # ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    # DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    # SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    # CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    # LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
    # OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
    # DAMAGE.

    def _m_values(self, j):
        size = 2 * j + 1
        if not size.is_integer() or not size > 0:
            raise ValueError(f"Only integer or half-integer values allowed for j, got: {j}")
        return size, [j - i for i in range(int(2 * j + 1))]

    def _eval_wignerd(self, j, m, mp, beta):
        from numpy import cos, pi, sin, sqrt

        try:
            from scipy.misc import factorial
        except ImportError:
            from scipy.special import factorial
        from scipy.special import binom as binomial

        r = 0
        if beta == pi / 2:
            # Varshalovich Equation (5), Section 4.16, page 113, setting
            # alpha=gamma=0.
            for k in range(int(2 * j) + 1):
                if k > j + mp or k > j - m or k < mp - m:
                    continue
                r += (-1) ** k * binomial(j + mp, k) * binomial(j - mp, k + m - mp)
            r *= (
                (-1) ** (m - mp)
                / 2**j
                * sqrt(factorial(j + m) * factorial(j - m) / (factorial(j + mp) * factorial(j - mp)))
            )
        else:
            # Varshalovich Equation(5), Section 4.7.2, page 87, where we set
            # beta1=beta2=pi/2, and we get alpha=gamma=pi/2 and beta=phi+pi,
            # then we use the Eq. (1), Section 4.4. page 79, to simplify:
            # d(j, m, mp, beta+pi) = (-1)**(j-mp) * d(j, m, -mp, beta)
            # This happens to be almost the same as in Eq.(10), Section 4.16,
            # except that we need to substitute -mp for mp.
            size, mvals = self._m_values(j)
            for mpp in mvals:
                r += (
                    self._eval_wignerd(j, m, mpp, pi / 2)
                    * (cos(-mpp * beta) + 1j * sin(-mpp * beta))
                    * self._eval_wignerd(j, mpp, -mp, pi / 2)
                )
            # Empirical normalization factor so results match Varshalovich
            # Tables 4.3-4.12
            # Note that this exact normalization does not follow from the
            # above equations
            r = r * 1j ** (2 * j - m - mp) * (-1) ** (2 * m)

        return r


# === Append csr matrices ===
# see http://stackoverflow.com/questions/4695337/expanding-adding-a-row-or-column-a-scipy-sparse-matrix


def csr_vappend(a, b) -> None:
    """Take in 2 csr_matrices and append the second one to the bottom of the first one.

    Much faster than scipy.sparse.vstack but assumes the type to be csr and overwrites
    the first matrix instead of copying it. The data, indices, and indptr still get copied.
    """
    if a.shape[1] != b.shape[1]:
        raise ValueError("Dimension mismatch in csr_vappend")

    a.data = np.hstack((a.data, b.data))
    a.indices = np.hstack((a.indices, b.indices))
    a.indptr = np.hstack((a.indptr, (b.indptr + a.nnz)[1:]))
    a._shape = (a.shape[0] + b.shape[0], b.shape[1])


# === Append csc matrices ===
# see http://stackoverflow.com/questions/4695337/expanding-adding-a-row-or-column-a-scipy-sparse-matrix


def csc_happend(a, b) -> None:
    """Take in 2 csc_matrices and append the second one to the right of the first one."""
    if a.shape[0] != b.shape[0]:
        raise ValueError("Dimension mismatch in csc_happend")

    a.data = np.hstack((a.data, b.data))
    a.indices = np.hstack((a.indices, b.indices))
    a.indptr = np.hstack((a.indptr, (b.indptr + a.nnz)[1:]))
    a._shape = (b.shape[0], a.shape[1] + b.shape[1])


# === Keep only the values that are maximal within a row of a csr matrix ===
# see http://stackoverflow.com/questions/15992857/efficient-way-to-get-the-max-of-each-row-for-large-sparse-matrix


def csr_keepmax(a) -> None:
    boolarr = np.diff(a.indptr) > 0
    ret = np.maximum.reduceat(a.data, a.indptr[:-1][boolarr])
    a.data[a.data != np.repeat(ret, np.diff(a.indptr)[boolarr])] = 0
    a.eliminate_zeros()


# === Return a normalized image ===


def normscale(data, cmin=None, cmax=None):
    if cmin is None:
        cmin = np.nanmin(data)
    if cmax is None:
        cmax = np.nanmax(data)
    return (data - cmin) / (cmax - cmin or 1)


# === Return a byte-scaled image ===


def bytescale(data, cmin=None, cmax=None, high=255, low=0):
    return np.array(normscale(data, cmin, cmax) * (high - low + 0.9999) + low).astype(int)
