# FieldCombinationTest
using Test
using PairInteraction

# test_combined_fields
cache = PairInteraction.MatrixElementCache()

# Set up SystemOne
system_one = PairInteraction.SystemOne("Rb", cache)
PairInteraction.restrictEnergy(system_one, -1077.243011609127, -939.9554235203701)
PairInteraction.restrictN(system_one, 57, 63)
PairInteraction.restrictL(system_one, 0, 3)

PairInteraction.setConservedMomentaUnderRotation(system_one, [-0.5f0])
PairInteraction.setEfield(system_one, [0, 0, 0.7])
PairInteraction.setBfield(system_one, [0, 0, -8.8])

# Diagonalize the system
PairInteraction.diagonalize(system_one)

# Compare results
hamiltonian = sparse(PairInteraction.getHamiltonian(system_one))
energies = [hamiltonian[i, i] for i in 1:size(hamiltonian)[1]]
@test isapprox(energies[14], -1000.26793709, atol=1e-4)
