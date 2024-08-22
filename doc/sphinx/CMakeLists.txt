# SPDX-FileCopyrightText: 2024 Sebastian Weber, Henri Menke, and contributors. All rights reserved.
# SPDX-License-Identifier: GPL-3.0-or-later OR LGPL-3.0-or-later

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

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${SPHINX_BASE}/_doctrees")

# HTML output directory
set(SPHINX_HTML_DIR "${SPHINX_BASE}/_build/html")

# Sphinx configuration file
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in" "${SPHINX_BASE}/conf.py" @ONLY)

# Directories to be copied to the build directory
file(COPY "_static" DESTINATION "${SPHINX_BASE}")
file(COPY "_templates" DESTINATION "${SPHINX_BASE}")

# Prepare sphinx target
add_custom_target(
  sphinx
  COMMAND # sphinx-build
          "${SPHINX_BUILD_EXECUTABLE}" -q -T -b html -c "${SPHINX_BASE}" -d "${SPHINX_CACHE_DIR}"
          "${CMAKE_CURRENT_SOURCE_DIR}" "${SPHINX_HTML_DIR}"
  VERBATIM
  DEPENDS bindings
  COMMENT "Building HTML documentation with Sphinx")

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples_python/ DESTINATION share/doc/pairinteraction/examples_python)
install(DIRECTORY ${SPHINX_HTML_DIR} DESTINATION share/doc/pairinteraction/html)

# Prepare livehtml target
add_custom_target(
  livehtml
  COMMAND # sphinx-autobuild
          "${SPHINX_AUTOBUILD_EXECUTABLE}" --ignore "*/_autosummary/*" -b html -c "${SPHINX_BASE}" -d
          "${SPHINX_CACHE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" "${SPHINX_HTML_DIR}"
  VERBATIM
  DEPENDS bindings
  COMMENT "Show a live preview of the documentation with Sphinx, serving on http://127.0.0.1:8000")

set_target_properties(livehtml PROPERTIES EXCLUDE_FROM_ALL ON)

# IPython notebooks to be tested
set(NOTEBOOKS "quick_start.ipynb")

foreach(file ${NOTEBOOKS})
  string(REGEX REPLACE "\\.[^.]*$" "" NOTEBOOK_NAME ${file})
  add_test(
    NAME "notebook.${NOTEBOOK_NAME}"
    COMMAND ${JUPYTER_NBCONVERT_EXECUTABLE} --to notebook --execute --ExecutePreprocessor.timeout=300
            "${CMAKE_CURRENT_SOURCE_DIR}/examples_python/${file}"
    WORKING_DIRECTORY "${SPHINX_BASE}")
  set_tests_properties(
    "notebook.${NOTEBOOK_NAME}"
    PROPERTIES
      ENVIRONMENT
      "OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;PYTHONPATH=$<SHELL_PATH:${CMAKE_BINARY_DIR};${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE};$ENV{PYTHONPATH}>"
  )
endforeach()
