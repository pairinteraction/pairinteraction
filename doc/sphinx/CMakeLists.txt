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
  find_package(IPython REQUIRED)

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
      ${PYTHON_EXECUTABLE} -c
      "import sys; from sphinx.apidoc import main; sys.exit(main(sys.argv))"
        -f
        -P
        -o ${SPHINX_BASE}
        ${SPHINX_PYTHON_DIR}
        ${EXCLUDE}
    COMMAND
      # override automatically generated modules.rst
      ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/modules.rst ${SPHINX_BASE}/modules.rst
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

  # Files to be copied to the build directory
  set(FILE_LIST
    "_templates/layout.html"
    "images/logo48x48.png"
    "html/slides/screen-osx-pairpotential.jpg"
    "html/slides/screen-osx-starkmap.jpg"
    "html/slides/screen-osx-zeemanmap.jpg"
    "html/slides/screen-win64.png"
    "html/buttons/github.png"
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
    "pair_potential_efield_sensitivity.ipynb")

  foreach(file ${NOTEBOOKS})
    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/examples_python/${file}"
      "${SPHINX_BASE}/${file}"
      COPYONLY)
    add_custom_target("${file}"
      COMMAND ${IPYTHON_EXECUTABLE} nbconvert --to rst "${SPHINX_BASE}/${file}"
      WORKING_DIRECTORY "${SPHINX_BASE}")
    add_dependencies(sphinx "${file}")
  endforeach()

endif(SPHINX_FOUND) 
