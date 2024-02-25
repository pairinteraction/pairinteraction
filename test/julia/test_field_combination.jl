# Copyright (c) 2020 Sebastian Weber, Henri Menke, Alexander Papageorge. All rights reserved.
#
# This file is part of the pairinteraction library.
#
# The pairinteraction library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.

# FieldCombinationTest
using Test
using PairInteraction

# test_combined_fields
cache = PairInteraction.MatrixElementCache()

# Set up SystemOne
system_one = PairInteraction.SystemOne{Float64}("Rb", cache)
PairInteraction.restrictEnergy(system_one, -1077.243011609127, -939.9554235203701)
PairInteraction.restrictN(system_one, 57, 63)
PairInteraction.restrictL(system_one, 0, 3)

PairInteraction.setConservedMomentaUnderRotation(system_one, [-0.5f0])
PairInteraction.setEfield(system_one, [0, 0, 0.7])
PairInteraction.setBfield(system_one, [0, 0, -8.8])
PairInteraction.enableDiamagnetism(system_one, false)

# Diagonalize the system
PairInteraction.diagonalize(system_one)

# Compare results
hamiltonian = sparse(PairInteraction.getHamiltonian(system_one))
energies = [hamiltonian[i, i] for i in 1:size(hamiltonian)[1]]
@test isapprox(energies[14], -1000.2679341660352, atol=1e-4)
