# Copyright (c) 2023 Sebastian Weber, Henri Menke. All rights reserved.
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

set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
set(BENCHMARK_ENABLE_TESTING OFF)

include(FetchContent)
FetchContent_Declare(
  GoogleBenchmark
  GIT_REPOSITORY https://github.com/google/benchmark
  GIT_TAG        v1.7.1
  GIT_PROGRESS   TRUE
)
FetchContent_MakeAvailable(GoogleBenchmark)

add_custom_target(run_benchmarks)

function(benchmark)
  cmake_parse_arguments(BENCH "" "TARGET" "SOURCE" ${ARGN})
  add_executable(${BENCH_TARGET} ${BENCH_SOURCE})
  target_compile_features(${BENCH_TARGET} PRIVATE cxx_std_17)
  set_target_properties(${BENCH_TARGET} PROPERTIES CXX_EXTENSIONS OFF)
  target_link_libraries(${BENCH_TARGET} PRIVATE pairinteraction benchmark::benchmark)
  add_custom_target(run-${BENCH_TARGET}
    ${BENCH_TARGET} --benchmark_format=json --benchmark_out=${BENCH_TARGET}.json
    DEPENDS ${BENCH_TARGET}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  add_dependencies(run_benchmarks run-${BENCH_TARGET})
endfunction()

benchmark(TARGET integration_bench SOURCE integration_bench.cpp)

if(WITH_LAPACK)
  benchmark(TARGET diagonalization_bench SOURCE diagonalization_bench.cpp)
  if(MKL_FOUND)
    target_link_libraries(diagonalization_bench PRIVATE MKL::MKL)
  elseif(LAPACK_FOUND)
    target_link_libraries(diagonalization_bench PRIVATE LAPACK::LAPACK LAPACKE::LAPACKE)
  endif()
endif()
