# Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

find_package(Sphinx)

if(SPHINX_FOUND)
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

  # Files to be copied to the build directory
  set(FILE_LIST 
    "index.rst")

  foreach(file ${FILE_LIST})
    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
      "${SPHINX_BASE}/${file}"
      COPYONLY)
  endforeach()

  add_custom_target(sphinx
      # sphinx-apidoc
      ${PYTHON_EXECUTABLE} -c
      "import sys; from sphinx.apidoc import main; sys.exit(main(sys.argv))"
        -f
        -P
        -o ${SPHINX_BASE}
        ${SPHINX_PYTHON_DIR}
        ${EXCLUDE}
    COMMAND
      # sphinx-build
      ${PYTHON_EXECUTABLE} -c
      "import sys; from sphinx import main; sys.exit(main(sys.argv))"
        -q
        -b html
        -c "${SPHINX_BASE}"
        -d "${SPHINX_CACHE_DIR}"
        "${CMAKE_CURRENT_BINARY_DIR}"
        "${SPHINX_HTML_DIR}"
    VERBATIM
    DEPENDS _pireal _picomplex
    COMMENT "Building HTML documentation with Sphinx")
endif(SPHINX_FOUND) 
