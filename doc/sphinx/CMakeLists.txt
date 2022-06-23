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
find_program(SPHINX_APIDOC_EXECUTABLE sphinx-apidoc REQUIRED)
if(NOT SPHINX_APIDOC_EXECUTABLE)
  message(FATAL_ERROR "sphinx-apidoc not found")
endif()
find_program(JUPYTER_NBCONVERT_EXECUTABLE jupyter-nbconvert REQUIRED)
if(NOT JUPYTER_NBCONVERT_EXECUTABLE)
  message(FATAL_ERROR "jupyter-nbconvert not found")
endif()

set(SPHINX_BASE "${CMAKE_CURRENT_BINARY_DIR}")

set(BINARY_BUILD_DIR "${SPHINX_BASE}/_build")
set(SPHINX_PYTHON_DIR "${CMAKE_BINARY_DIR}/${LIBNAME}")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${SPHINX_BASE}/_doctrees")

# HTML output directory
set(SPHINX_HTML_DIR "${SPHINX_BASE}/html")

# Sphinx configuration file
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
  "${SPHINX_BASE}/conf.py"
  @ONLY)

# Prepare sphinx target
add_custom_target(sphinx
    ${CMAKE_COMMAND} -E env SPHINX_APIDOC_OPTIONS=members,inherited-members,show-inheritance
    # sphinx-apidoc
    "${SPHINX_APIDOC_EXECUTABLE}"
      -f
      -P
      -o ${SPHINX_BASE}
      ${SPHINX_PYTHON_DIR}
      ${EXCLUDE}
  COMMAND
    # override automatically generated modules.rst
    ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/modules.rst ${SPHINX_BASE}/modules.rst
  COMMAND
    # override automatically generated pairinteraction.rst
    ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pairinteraction.rst ${SPHINX_BASE}/pairinteraction.rst
  COMMAND
    # sphinx-build
    "${SPHINX_BUILD_EXECUTABLE}"
      -q
      -b html
      -c "${SPHINX_BASE}"
      -d "${SPHINX_CACHE_DIR}"
      "${CMAKE_CURRENT_BINARY_DIR}"
      "${SPHINX_HTML_DIR}"
  VERBATIM
  DEPENDS _pireal _picomplex
  COMMENT "Building HTML documentation with Sphinx")

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples_cpp/
  DESTINATION share/doc/pairinteraction/examples_cpp)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples_python/
  DESTINATION share/doc/pairinteraction/examples_python)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
  DESTINATION share/doc/pairinteraction/html)

# Files to be copied to the build directory
set(FILE_LIST
  "_templates/layout.html"
  "images/logo48x48.png"
  "html/slides/screen-osx-pairpotential.jpg"
  "html/slides/screen-osx-starkmap.jpg"
  "html/slides/screen-osx-zeemanmap.jpg"
  "html/slides/screen-win64.png"
  "html/buttons/github.svg"
  "images/license-badge.svg"
  "images/arXiv-badge.svg"
  "index.rst"
  "installation.rst"
  "genindex.rst"
  "tutorials.rst")

foreach(file ${FILE_LIST})
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
    "${SPHINX_BASE}/${file}"
    COPYONLY)
endforeach()

# IPython notebooks to be processed
set(NOTEBOOKS
  "introduction.ipynb"
  "matrix_elements.ipynb"
  "wavefunctions.ipynb"
  "comparison_to_saffman_fig13.ipynb"
  "pair_potential_efield_sensitivity.ipynb"
  "pair_potential_near_surface.ipynb"
  "vdw_near_surface.ipynb"
  "atom_ion_interaction.ipynb")

foreach(file ${NOTEBOOKS})
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/examples_python/${file}"
    "${SPHINX_BASE}/${file}"
    COPYONLY)
  add_custom_target("${file}"
    COMMAND ${JUPYTER_NBCONVERT_EXECUTABLE} --to rst "${SPHINX_BASE}/${file}"
    WORKING_DIRECTORY "${SPHINX_BASE}")
  add_dependencies(sphinx "${file}")
endforeach()

# IPython notebooks to be tested
set(NOTEBOOKS
  "introduction.ipynb"
  "matrix_elements.ipynb"
  "vdw_near_surface.ipynb"
  "atom_ion_interaction.ipynb")

if(WITH_GSL)
  list(APPEND NOTEBOOKS "wavefunctions.ipynb")
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  list(APPEND NOTEBOOKS
    "comparison_to_saffman_fig13.ipynb"
    "pair_potential_efield_sensitivity.ipynb"
    "pair_potential_near_surface.ipynb")
endif()

foreach(file ${NOTEBOOKS})
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/examples_python/${file}"
    "${SPHINX_BASE}/${file}"
    COPYONLY)
  string(REGEX REPLACE "\\.[^.]*$" "" NOTEBOOK_NAME ${file})
  add_test(NAME "notebook.${NOTEBOOK_NAME}"
    COMMAND ${JUPYTER_NBCONVERT_EXECUTABLE} --execute --ExecutePreprocessor.timeout=300 "${SPHINX_BASE}/${file}"
    WORKING_DIRECTORY "${SPHINX_BASE}")
  set_tests_properties("notebook.${NOTEBOOK_NAME}"
    PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR};${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}")
endforeach()
