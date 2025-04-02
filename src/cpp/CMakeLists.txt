# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

add_library(
  pairinteraction OBJECT
  # Headers
  ./include/pairinteraction/basis/Basis.hpp
  ./include/pairinteraction/basis/BasisAtom.hpp
  ./include/pairinteraction/basis/BasisAtomCreator.hpp
  ./include/pairinteraction/basis/BasisPair.hpp
  ./include/pairinteraction/basis/BasisPairCreator.hpp
  ./include/pairinteraction/database/AtomDescriptionByParameters.hpp
  ./include/pairinteraction/database/AtomDescriptionByRanges.hpp
  ./include/pairinteraction/database/Database.hpp
  ./include/pairinteraction/database/GitHubDownloader.hpp
  ./include/pairinteraction/database/ParquetManager.hpp
  ./include/pairinteraction/diagonalize/diagonalize.hpp
  ./include/pairinteraction/diagonalize/DiagonalizerEigen.hpp
  ./include/pairinteraction/diagonalize/DiagonalizerFeast.hpp
  ./include/pairinteraction/diagonalize/DiagonalizerLapackeEvd.hpp
  ./include/pairinteraction/diagonalize/DiagonalizerLapackeEvr.hpp
  ./include/pairinteraction/enums/FloatType.hpp
  ./include/pairinteraction/enums/OperatorType.hpp
  ./include/pairinteraction/enums/Parity.hpp
  ./include/pairinteraction/enums/TransformationType.hpp
  ./include/pairinteraction/interfaces/DiagonalizerInterface.hpp
  ./include/pairinteraction/interfaces/TransformationBuilderInterface.hpp
  ./include/pairinteraction/ket/Ket.hpp
  ./include/pairinteraction/ket/KetAtom.hpp
  ./include/pairinteraction/ket/KetAtomCreator.hpp
  ./include/pairinteraction/ket/KetPair.hpp
  ./include/pairinteraction/operator/Operator.hpp
  ./include/pairinteraction/operator/OperatorAtom.hpp
  ./include/pairinteraction/operator/OperatorPair.hpp
  ./include/pairinteraction/system/GreenTensor.hpp
  ./include/pairinteraction/system/System.hpp
  ./include/pairinteraction/system/SystemAtom.hpp
  ./include/pairinteraction/system/SystemPair.hpp
  ./include/pairinteraction/tools/run_unit_tests.hpp
  ./include/pairinteraction/tools/setup.hpp
  ./include/pairinteraction/utils/args.hpp
  ./include/pairinteraction/utils/eigen_assertion.hpp
  ./include/pairinteraction/utils/eigen_compat.hpp
  ./include/pairinteraction/utils/euler.hpp
  ./include/pairinteraction/utils/hash.hpp
  ./include/pairinteraction/utils/id_in_database.hpp
  ./include/pairinteraction/utils/maths.hpp
  ./include/pairinteraction/utils/paths.hpp
  ./include/pairinteraction/utils/Range.hpp
  ./include/pairinteraction/utils/spherical.hpp
  ./include/pairinteraction/utils/streamed.hpp
  ./include/pairinteraction/utils/tensor.hpp
  ./include/pairinteraction/utils/traits.hpp
  ./include/pairinteraction/utils/wigner.hpp
  ./include/pairinteraction/pairinteraction.hpp
  ./include/pairinteraction/Info.hpp
  ./include/pairinteraction/version.hpp
  # Sources
  ./src/basis/Basis.cpp
  ./src/basis/BasisAtom.cpp
  ./src/basis/BasisAtomCreator.cpp
  ./src/basis/BasisAtomCreator.test.cpp
  ./src/basis/BasisPair.cpp
  ./src/basis/BasisPairCreator.cpp
  ./src/basis/BasisPairCreator.test.cpp
  ./src/database/Database.cpp
  ./src/database/Database.test.cpp
  ./src/database/GitHubDownloader.cpp
  ./src/database/GitHubDownloader.test.cpp
  ./src/database/ParquetManager.cpp
  ./src/database/ParquetManager.test.cpp
  ./src/diagonalize/diagonalize.cpp
  ./src/diagonalize/DiagonalizerEigen.cpp
  ./src/diagonalize/DiagonalizerFeast.cpp
  ./src/diagonalize/DiagonalizerLapackeEvd.cpp
  ./src/diagonalize/DiagonalizerLapackeEvr.cpp
  ./src/enums/Parity.test.cpp
  ./src/interfaces/DiagonalizerInterface.cpp
  ./src/interfaces/TransformationBuilderInterface.cpp
  ./src/ket/Ket.cpp
  ./src/ket/Ket.test.cpp
  ./src/ket/KetAtom.cpp
  ./src/ket/KetAtomCreator.cpp
  ./src/ket/KetAtomCreator.test.cpp
  ./src/ket/KetPair.cpp
  ./src/operator/Operator.cpp
  ./src/operator/OperatorAtom.cpp
  ./src/operator/OperatorPair.cpp
  ./src/system/GreenTensor.cpp
  ./src/system/System.cpp
  ./src/system/SystemAtom.cpp
  ./src/system/SystemAtom.test.cpp
  ./src/system/SystemPair.cpp
  ./src/system/SystemPair.test.cpp
  ./src/tools/setup.cpp
  ./src/tools/run_unit_tests.cpp
  ./src/utils/euler.test.cpp
  ./src/utils/spherical.cpp
  ./src/utils/spherical.test.cpp
  ./src/utils/tensor.cpp
  ./src/utils/wigner.test.cpp)

target_link_libraries(
  pairinteraction
  PUBLIC $<$<BOOL:${WITH_COVERAGE}>:coverage::gcov>
         doctest::doctest
         duckdb
         Eigen3::Eigen
         httplib::httplib
         nlohmann_json::nlohmann_json
         spdlog::spdlog
         TBB::tbb
         fmt::fmt
         cpptrace::cpptrace
         miniz
         $<$<TARGET_EXISTS:MKL::MKL>:MKL::MKL>
         $<$<AND:$<TARGET_EXISTS:LAPACKE::LAPACKE>,$<NOT:$<TARGET_EXISTS:MKL::MKL>>>:LAPACKE::LAPACKE>)

target_include_directories(pairinteraction PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_compile_definitions(
  pairinteraction
  PUBLIC DOCTEST_CONFIG_SUPER_FAST_ASSERTS
         $<$<TARGET_EXISTS:MKL::MKL>:WITH_MKL
         EIGEN_USE_MKL_ALL
         EIGEN_MKL_NO_DIRECT_CALL
         WITH_LAPACKE>
         $<$<AND:$<TARGET_EXISTS:LAPACKE::LAPACKE>,$<NOT:$<TARGET_EXISTS:MKL::MKL>>>:WITH_LAPACKE
         EIGEN_USE_BLAS
         EIGEN_USE_LAPACKE>
         SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG
         $<$<PLATFORM_ID:Windows>:WIN32_LEAN_AND_MEAN>)

if(TARGET MKL::MKL)
  target_compile_definitions(pairinteraction PUBLIC MKL_Complex8=std::complex<float> MKL_Complex16=std::complex<double>)
elseif(TARGET LAPACKE::LAPACKE)
  target_compile_definitions(pairinteraction PUBLIC lapack_complex_float=std::complex<float>
                                                    lapack_complex_double=std::complex<double>)
endif()

target_compile_options(
  pairinteraction
  PUBLIC $<$<CXX_COMPILER_ID:MSVC>:/W4 /bigobj /MP>
         $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra
         -Wpedantic>)
target_compile_features(pairinteraction PUBLIC cxx_std_17)
