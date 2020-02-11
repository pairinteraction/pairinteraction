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

# TestPairState
using Test
using PairInteraction

# setUp
s1 = PairInteraction.StateOne("Sr3", 80, 2, 1f0, 0f0)
s2 = PairInteraction.StateOne("Sr3", 79, 1, 2f0, 0f0)
s = PairInteraction.StateTwo(s1, s2)
s_artificial = PairInteraction.StateTwo(["G", "G"])

# test_combination
@test typeof(s) === PairInteraction.StateTwoAllocated

# test_properties
@test PairInteraction.getFirstState(s) == s1
@test PairInteraction.getSecondState(s) == s2
@test PairInteraction.getSpecies(s) == ("Sr3", "Sr3")
@test PairInteraction.getElement(s) == ("Sr", "Sr")
@test PairInteraction.getS(s) == (1, 1)
@test PairInteraction.getN(s) == (80, 79)
@test PairInteraction.getL(s) == (2, 1)
@test PairInteraction.getJ(s) == (1, 2)
@test PairInteraction.getM(s) == (0, 0)
@test PairInteraction.getLabel(s_artificial) == ("G", "G")
@test PairInteraction.isArtificial(s_artificial) == (true, true)

for i in (0, 1)
    @test PairInteraction.getSpecies(s, i) == ("Sr3", "Sr3")[i+1]
    @test PairInteraction.getElement(s, i) == ("Sr", "Sr")[i+1]
    @test PairInteraction.getS(s, i) == [1, 1][i+1]
    @test PairInteraction.getN(s, i) == [80, 79][i+1]
    @test PairInteraction.getL(s, i) == [2, 1][i+1]
    @test PairInteraction.getJ(s, i) == [1, 2][i+1]
    @test PairInteraction.getM(s, i) == [0, 0][i+1]
    @test PairInteraction.getLabel(s_artificial, i) == ("G", "G")[i+1]
    @test PairInteraction.isArtificial(s_artificial, i) == (true, true)[i+1]
end
