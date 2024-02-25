#include "MatrixElementCache.hpp"
#include "State.hpp"
#include "SystemOne.hpp"
#include "SystemTwo.hpp"
#include "filesystem.hpp"

#include <benchmark/benchmark.h>

static void benchmark_integration(benchmark::State &state) {
    using Scalar = double;

    // Set up cache
    fs::path path_cache = fs::create_temp_directory();
    MatrixElementCache cache(path_cache);

    for (auto _ : state) {
        // Setup states
        StateOne state_one("Rb", 61, 2, 1.5, 1.5);
        StateTwo state_two(state_one, state_one);

        ////////////////////////////////////////////////////////////////////
        /// Test system consisting of one atom /////////////////////////////
        ////////////////////////////////////////////////////////////////////

        // Build one-atom system
        SystemOne<Scalar> system_one(state_one.getSpecies(), cache);
        system_one.restrictEnergy(state_one.getEnergy() - 40, state_one.getEnergy() + 40);
        system_one.restrictN(state_one.getN() - 1, state_one.getN() + 1);
        system_one.restrictL(state_one.getL() - 1, state_one.getL() + 1);
        system_one.setEfield({{0, 0, 0.1}});
        system_one.setBfield({{0, 0, 1}});
        system_one.enableDiamagnetism(false);

        // Diagonalize one-atom system
        system_one.diagonalize();

        ////////////////////////////////////////////////////////////////////
        /// Test system consisting of two atoms ////////////////////////////
        ////////////////////////////////////////////////////////////////////

        // Build one-atom system (for this test, system_one has to be diagonal by
        // itself because diagonalization can lead to different order of
        // eigenvectors)
        // system_one = SystemOne(state_one.species, cache); // TODO  object of type 'SystemOne'
        // cannot be assigned because its copy assignment operator is implicitly deleted
        SystemOne<Scalar> system_one_new(state_one.getSpecies(), cache);
        system_one_new.restrictEnergy(state_one.getEnergy() - 40, state_one.getEnergy() + 40);
        system_one_new.restrictN(state_one.getN() - 1, state_one.getN() + 1);
        system_one_new.restrictL(state_one.getL() - 1, state_one.getL() + 1);

        // Build two-atom system
        SystemTwo<Scalar> system_two(system_one_new, system_one_new, cache);
        system_two.restrictEnergy(state_two.getEnergy() - 2, state_two.getEnergy() + 2);
        system_two.setConservedParityUnderPermutation(ODD);
        system_two.setDistance(6);
        system_two.setAngle(0.9);

        // Diagonalize two-atom system
        system_two.diagonalize();
    }

    // Delete cache directory
    fs::remove_all(path_cache);
}

BENCHMARK(benchmark_integration);

BENCHMARK_MAIN();
