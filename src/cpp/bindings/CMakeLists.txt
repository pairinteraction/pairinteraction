# Add bindings object
add_library(bindings OBJECT Ket.py.cpp Ket.py.hpp binding.cpp)

target_include_directories(
  bindings
  PRIVATE "$<TARGET_PROPERTY:pairinteraction,INCLUDE_DIRECTORIES>;${Python_INCLUDE_DIRS};${Python_NumPy_INCLUDE_DIRS}")
target_link_libraries(bindings PRIVATE $<$<BOOL:${WITH_COVERAGE}>:coverage::gcov> nanobind-static)
target_compile_options(bindings PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/bigobj /MP>)

# Assemble objects into backend library
add_library(backend SHARED)

target_link_libraries(
  backend
  PUBLIC $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARIES_GENERALIZED}>
  PRIVATE pairinteraction bindings)

set_target_properties(
  backend
  PROPERTIES PREFIX ""
             SUFFIX "${NB_SUFFIX_EXT}"
             INSTALL_RPATH $<$<BOOL:${SKBUILD}>:$<$<PLATFORM_ID:Linux>:$ORIGIN>$<$<PLATFORM_ID:Darwin>:@loader_path>>)

target_link_options(
  backend
  PRIVATE
  $<$<CXX_COMPILER_ID:AppleClang>:-undefined
  dynamic_lookup>
  $<$<CXX_COMPILER_ID:Clang,GNU>:-Wl,--unresolved-symbols=ignore-all>
  $<$<CXX_COMPILER_ID:MSVC>:/NODEFAULTLIB:python3${Python_VERSION_MINOR}.lib
  /DEFAULTLIB:python3.lib
  /FORCE:UNRESOLVED>)

if(WIN32)
  add_custom_command(
    TARGET backend
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:backend> $<TARGET_FILE_DIR:backend>
    COMMAND_EXPAND_LISTS)
endif()

install(IMPORTED_RUNTIME_ARTIFACTS duckdb DESTINATION pairinteraction)
install(TARGETS backend DESTINATION pairinteraction)
if(WIN32)
  install(FILES $<TARGET_RUNTIME_DLLS:backend> DESTINATION pairinteraction)
endif()
