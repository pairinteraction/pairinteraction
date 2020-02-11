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

# TestState
using Test
using PairInteraction

# setUp
s = PairInteraction.StateOne("Sr3", 79, 1, 2.f0, 0.f0)
s_artificial = PairInteraction.StateOne("G")

# test_instantiation
@test typeof(s) === PairInteraction.StateOneAllocated

# test_properties
@test PairInteraction.getSpecies(s) == "Sr3"
@test PairInteraction.getElement(s) == "Sr"
@test PairInteraction.getS(s) == 1
@test PairInteraction.getN(s) == 79
@test PairInteraction.getL(s) == 1
@test PairInteraction.getJ(s) == 2
@test PairInteraction.getM(s) == 0
@test PairInteraction.getLabel(s_artificial) == "G"
@test PairInteraction.isArtificial(s_artificial) === true

# test_comparison
@test s == PairInteraction.StateOne("Sr3", 79, 1, 2.f0, 0.f0)
@test (s != PairInteraction.StateOne("Sr3", 79, 1, 2.f0, 0.f0)) == false
@test s != PairInteraction.StateOne("Sr3", 79, 1, 2.f0, 1.f0)
@test (s == PairInteraction.StateOne("Sr3", 79, 1, 2.f0, 1.f0)) == false
@test s ^ PairInteraction.StateOne("Sr3", 79, 1, 2.f0, 0.f0)
@test (s ^ PairInteraction.StateOne("Sr3", 79, 2, 2.f0, 0.f0)) == false
@test s ^ PairInteraction.StateOne("Sr3", PairInteraction.ARB, 1, 2.f0, 0.f0)
@test (s ^ PairInteraction.StateOne("Sr3", PairInteraction.ARB, 2, 2.f0, 0.f0)) == false
@test (s == s_artificial) == false

# test_output
@test string(s) == "|Sr3, 79 P_2, mj=0>"

# test_reflection
s0 = PairInteraction.StateOne("Sr3", 79, 1, 2.f0, 1.f0)
s1 = PairInteraction.StateOne("Sr3", 79, 1, 2.f0, -1.f0)
@test PairInteraction.getReflected(s0) == s1
