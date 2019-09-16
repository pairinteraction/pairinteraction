# GreenTensorTest

using Test
using SparseArrays: dropzeros!
using PairInteraction

# setUp
# Set up cache
cache = PairInteraction.MatrixElementCache()

# Setup states
state_one = PairInteraction.StateOne("Rb", 61, 1, 1.5f0, 1.5f0)
state_two = PairInteraction.StateTwo(state_one, state_one)

# test_greentensor_angle
# Build one-atom system
system_one = PairInteraction.SystemOne(PairInteraction.getSpecies(state_one), cache)
PairInteraction.restrictEnergy(system_one, PairInteraction.getEnergy(state_one) - 40, 
                                           PairInteraction.getEnergy(state_one) + 40)
PairInteraction.restrictN(system_one, PairInteraction.getN(state_one) - 1, 
                                      PairInteraction.getN(state_one) + 1)
PairInteraction.restrictL(system_one, PairInteraction.getL(state_one) - 2, 
                                      PairInteraction.getL(state_one) + 2)

# Build two-atom system
system_two = PairInteraction.SystemTwo(system_one, system_one, cache)
PairInteraction.restrictEnergy(system_two, PairInteraction.getEnergy(state_two) - 5, 
                                           PairInteraction.getEnergy(state_two) + 5)
PairInteraction.setConservedParityUnderInversion(system_two, PairInteraction.ODD)
PairInteraction.setConservedParityUnderPermutation(system_two, PairInteraction.ODD)
PairInteraction.setDistance(system_two, 5)
PairInteraction.setAngle(system_two, 1.78)

# Construct the Hamiltonian using the standard approach
system_two_standard = PairInteraction.SystemTwo(system_two)
PairInteraction.enableGreenTensor(system_two_standard, false)
hamiltonian_standard = sparse(PairInteraction.getHamiltonian(system_two_standard))

# Construct the Hamiltonian using the green tensor approach
system_two_greentensor = PairInteraction.SystemTwo(system_two)
PairInteraction.enableGreenTensor(system_two_greentensor, true)
hamiltonian_greentensor = sparse(PairInteraction.getHamiltonian(system_two_greentensor))

# Prune Hamiltonians (without pruning, max_diff_hamiltonian might be infinity due to division by zero)
hamiltonian_standard[abs.(hamiltonian_standard) .< 1e-6] .= 0
dropzeros!(hamiltonian_standard)
hamiltonian_greentensor[abs.(hamiltonian_greentensor) .< 1e-6] .= 0
dropzeros!(hamiltonian_greentensor)

# Compare Hamiltonians
@test isapprox(hamiltonian_standard, hamiltonian_greentensor)


# test_greentensor_surface
theta = π / 2
interatomic_distance = 10
distance_to_surface = [2.65 / 6, 5.29 / 6, 7.9 / 6] * interatomic_distance  # center of mass distance
state_one1 = PairInteraction.StateOne("Rb", 69, 0, 0.5f0, 0.5f0)
state_one2 = PairInteraction.StateOne("Rb", 72, 0, 0.5f0, 0.5f0)

# Set up pair state
state_two = PairInteraction.StateTwo(state_one1, state_one2)

# Set up one-atom system
system_one = PairInteraction.SystemOne(PairInteraction.getSpecies(state_one1), cache)
PairInteraction.restrictEnergy(system_one, 
                               min(PairInteraction.getEnergy(state_one1), 
                               PairInteraction.getEnergy(state_one2)) - 30,
                               max(PairInteraction.getEnergy(state_one1), 
                                   PairInteraction.getEnergy(state_one2)) + 30)
PairInteraction.restrictN(system_one, 
                          min(PairInteraction.getN(state_one1), PairInteraction.getN(state_one2)) - 2,
                          max(PairInteraction.getN(state_one1), PairInteraction.getN(state_one2)) + 2)
PairInteraction.restrictL(system_one,
                          min(PairInteraction.getL(state_one1), PairInteraction.getL(state_one2)) - 1,
                          max(PairInteraction.getL(state_one1), PairInteraction.getL(state_one2)) + 1)

# Set up two-atom system
system_two = PairInteraction.SystemTwo(system_one, system_one, cache)
PairInteraction.restrictEnergy(system_two, 
                               PairInteraction.getEnergy(state_two) - 3, 
                               PairInteraction.getEnergy(state_two) + 3)

PairInteraction.setAngle(system_two, theta)
PairInteraction.setDistance(system_two, interatomic_distance)
PairInteraction.enableGreenTensor(system_two, true)

# Calculate dispersion coefficients
PairInteraction.diagonalize(system_two)
idx = argmax(PairInteraction.getOverlap(system_two, state_two, 0, -theta, 0))

C6_freespace = (sparse(PairInteraction.getHamiltonian(system_two))[idx, idx] 
                - PairInteraction.getEnergy(state_two)) * interatomic_distance^6

PairInteraction.setSurfaceDistance(system_two, distance_to_surface[1])
PairInteraction.diagonalize(system_two)
idx = argmax(PairInteraction.getOverlap(system_two, state_two, 0, -theta, 0))
C6_1 = (sparse(PairInteraction.getHamiltonian(system_two))[idx, idx] 
        - PairInteraction.getEnergy(state_two)) * interatomic_distance^6

PairInteraction.setSurfaceDistance(system_two, distance_to_surface[2])
PairInteraction.diagonalize(system_two)
idx = argmax(PairInteraction.getOverlap(system_two, state_two, 0, -theta, 0))
C6_2 = (sparse(PairInteraction.getHamiltonian(system_two))[idx, idx] 
        - PairInteraction.getEnergy(state_two)) * interatomic_distance^6

PairInteraction.setSurfaceDistance(system_two, distance_to_surface[3])
PairInteraction.diagonalize(system_two)
idx = argmax(PairInteraction.getOverlap(system_two, state_two, 0, -theta, 0))
C6_3 = (sparse(PairInteraction.getHamiltonian(system_two))[idx, idx] 
        - PairInteraction.getEnergy(state_two)) * interatomic_distance^6

# Compare the results against previously calculated values
@test isapprox(real(C6_freespace), -670, atol=20)
@test isapprox(real(C6_1), -544, atol=20)
@test isapprox(real(C6_2), -628, atol=20)
@test isapprox(real(C6_3), -649, atol=20)
