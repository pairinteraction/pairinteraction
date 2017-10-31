include(FindPackageHandleStandardArgs)

if(NOT PYTHON_EXECUTABLE)
  find_package(PythonInterp 3)
endif()

if(PYTHON_EXECUTABLE)
  # find sphinx-build
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
    "import sphinx; print('TRUE')"
    OUTPUT_VARIABLE SPHINX_BUILD_FOUND
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # find sphinx-apidoc
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
    "import sphinx.apidoc; print('TRUE')"
    OUTPUT_VARIABLE SPHINX_API_DOC_FOUND
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Get the version
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
    "import sphinx; print(sphinx.__version__)"
    OUTPUT_VARIABLE SPHINX_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)

else()
  message(STATUS "Python not found, not finding Sphinx")
endif()

find_package_handle_standard_args(Sphinx 
    REQUIRED_VARS SPHINX_BUILD_FOUND SPHINX_API_DOC_FOUND
    VERSION_VAR SPHINX_VERSION
)
 
mark_as_advanced(SPHINX_EXECUTABLE)
