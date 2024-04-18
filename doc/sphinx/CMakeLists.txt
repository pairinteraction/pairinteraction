# Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
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

find_program(SPHINX_BUILD_EXECUTABLE sphinx-build REQUIRED)
if(NOT SPHINX_BUILD_EXECUTABLE)
  message(FATAL_ERROR "sphinx-build not found")
endif()
find_program(SPHINX_AUTOBUILD_EXECUTABLE sphinx-autobuild)
if(NOT SPHINX_AUTOBUILD_EXECUTABLE)
  message(STATUS "sphinx-autobuild not found")
else()
  message(STATUS "sphinx-autobuild found: ${SPHINX_AUTOBUILD_EXECUTABLE}")
endif()
find_program(JUPYTER_NBCONVERT_EXECUTABLE jupyter-nbconvert REQUIRED)
if(NOT JUPYTER_NBCONVERT_EXECUTABLE)
  message(FATAL_ERROR "jupyter-nbconvert not found")
endif()

set(SPHINX_BASE "${CMAKE_CURRENT_BINARY_DIR}")
set(BINARY_BUILD_DIR "${SPHINX_BASE}/_build")
set(SPHINX_PYTHON_DIR "${CMAKE_BINARY_DIR}/pairinteraction_backend_old")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${SPHINX_BASE}/_doctrees")

# HTML output directory
set(SPHINX_HTML_DIR "${SPHINX_BASE}/_build/html")

# Sphinx configuration file
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
  "${SPHINX_BASE}/conf.py"
  @ONLY)

# Directories to be copied to the build directory
file(COPY "_static" DESTINATION "${SPHINX_BASE}")
file(COPY "_templates" DESTINATION "${SPHINX_BASE}")

# Prepare sphinx target
add_custom_target(sphinx
  COMMAND
    # sphinx-build
    "${SPHINX_BUILD_EXECUTABLE}"
    -q
    -T
    -b html
    -c "${SPHINX_BASE}"
    -d "${SPHINX_CACHE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${SPHINX_HTML_DIR}"
  VERBATIM
  DEPENDS binding
  COMMENT "Building HTML documentation with Sphinx")

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples_cpp/
  DESTINATION share/doc/pairinteraction/examples_cpp)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples_python/
  DESTINATION share/doc/pairinteraction/examples_python)
install(DIRECTORY ${SPHINX_HTML_DIR}
  DESTINATION share/doc/pairinteraction/html)

# Prepare livehtml target
add_custom_target(livehtml
  COMMAND
    # sphinx-autobuild
    "${SPHINX_AUTOBUILD_EXECUTABLE}"
    --ignore "_autosummary"
    -b html
    -c "${SPHINX_BASE}"
    -d "${SPHINX_CACHE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${SPHINX_HTML_DIR}"
  VERBATIM
  DEPENDS binding
  COMMENT "Show a live preview of the documentation with Sphinx, serving on http://127.0.0.1:8000")

set_target_properties(livehtml PROPERTIES EXCLUDE_FROM_ALL ON)

# IPython notebooks to be tested
set(NOTEBOOKS
  "introduction.ipynb"
  "matrix_elements.ipynb"
  "vdw_near_surface.ipynb")

if(WITH_GSL)
  list(APPEND NOTEBOOKS "wavefunctions.ipynb")
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  list(APPEND NOTEBOOKS
    "comparison_to_saffman_fig13.ipynb"
    "pair_potential_efield_sensitivity.ipynb"
    "pair_potential_near_surface.ipynb"
    "atom_ion_interaction.ipynb")
endif()

foreach(file ${NOTEBOOKS})
  string(REGEX REPLACE "\\.[^.]*$" "" NOTEBOOK_NAME ${file})
  add_test(NAME "notebook.${NOTEBOOK_NAME}"
    COMMAND ${JUPYTER_NBCONVERT_EXECUTABLE} --to notebook --execute --ExecutePreprocessor.timeout=300 "${CMAKE_CURRENT_SOURCE_DIR}/examples_python/${file}"
    WORKING_DIRECTORY "${SPHINX_BASE}")
  set_tests_properties("notebook.${NOTEBOOK_NAME}"
    PROPERTIES ENVIRONMENT "OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;PYTHONPATH=$<SHELL_PATH:${CMAKE_BINARY_DIR};${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE};$ENV{PYTHONPATH}>"
)
endforeach()
