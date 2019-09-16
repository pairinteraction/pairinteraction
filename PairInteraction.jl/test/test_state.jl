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
