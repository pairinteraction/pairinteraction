# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
#
# This file is part of the pairinteraction GUI.
#
# The pairinteraction GUI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction GUI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the pairinteraction GUI. If not, see <http://www.gnu.org/licenses/>.
import locale

from pint import UnitRegistry
from scipy import constants

# === Initialize unit registry ===

locale.setlocale(locale.LC_NUMERIC, "C")

ureg = UnitRegistry()
Q = ureg.Quantity


def U(u):
    return ureg.Quantity(1, u)


def C(s):
    return Q(constants.value(s), constants.unit(s))


# === Define units that are used in the python program for the gui (the c++ program uses atomic units) ===


class Units:
    length = "micrometer"
    energy = "gigahertz"
    efield = "volt/centimeter"
    bfield = "gauss"
    angle = "degree"
    au_length = "au_length"
    au_energy = "au_energy"
    au_efield = "au_efield"
    au_bfield = "au_bfield"
    au_angle = "au_angle"


# === Handle quantities ===


class Quantity:
    au = {}
    au[Units.au_length] = C("atomic unit of length")
    au[Units.au_energy] = C("atomic unit of energy") / C("Planck constant")
    au[Units.au_efield] = C("atomic unit of electric field")
    au[Units.au_bfield] = C("atomic unit of mag. flux density")
    au[Units.au_angle] = Q("1")

    converter_au = {}
    converter_au[str(U(Units.length).dimensionality)] = [Units.au_length, au[Units.au_length]]
    converter_au[str(U(Units.energy).dimensionality)] = [Units.au_energy, au[Units.au_energy]]
    converter_au[str(U(Units.efield).dimensionality)] = [Units.au_efield, au[Units.au_efield]]
    converter_au[str(U(Units.bfield).dimensionality)] = [Units.au_bfield, au[Units.au_bfield]]
    converter_au[str(U(Units.angle).dimensionality)] = [Units.au_angle, au[Units.au_angle]]

    converter_uu = {}
    converter_uu[str(U(Units.length).dimensionality)] = Units.length
    converter_uu[str(U(Units.energy).dimensionality)] = Units.energy
    converter_uu[str(U(Units.efield).dimensionality)] = Units.efield
    converter_uu[str(U(Units.bfield).dimensionality)] = Units.bfield
    converter_uu[str(U(Units.angle).dimensionality)] = Units.angle

    def __init__(self, magnitude, units=None):
        self._magnitude = magnitude
        if units is None:
            self._units = units
        else:
            self._units = str(units)

    @property
    def magnitude(self):
        return self._magnitude

    @property
    def units(self):
        return self._units

    def toAU(self):
        if self.units is None:
            return self
            # raise Exception("Quantity can not be converted to atomic units.")

        if self.units[:3] == "au_":  # already in atomic units
            return self

        else:  # from other units to atomic units
            oldunits = U(self.units)
            newunits, converter = self.converter_au[str(oldunits.dimensionality)]

            if self.magnitude is None:
                return Quantity(self.magnitude, newunits)

            quantity = Q(self.magnitude, oldunits)

            quantity = (quantity / converter).to("dimensionless")
            return Quantity(quantity.magnitude, newunits)

    def toUU(self):
        if self.units is None:
            return self
            # raise Exception("Quantity can not be converted to user units.")

        if self.units[:3] == "au_":  # from atomic units to user units
            oldunits = self.au[self.units]
            newunits = self.converter_uu[str(oldunits.dimensionality)]

            if self.magnitude is None:
                return Quantity(self.magnitude, newunits)

            quantity = self.magnitude * oldunits

        else:  # from other units to user units
            oldunits = U(self.units)
            newunits = self.converter_uu[str(oldunits.dimensionality)]

            if self.magnitude is None:
                return Quantity(self.magnitude, newunits)

            quantity = Q(self.magnitude, oldunits)

        quantity = quantity.to(newunits)
        return Quantity(quantity.magnitude, newunits)
