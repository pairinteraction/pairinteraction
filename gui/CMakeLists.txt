# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
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

find_package(PyUIC REQUIRED)

configure_file(version.py.cmakein ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/version.py)

add_custom_target(gui ALL)

add_custom_command(
  TARGET gui
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/conf ${CMAKE_CURRENT_BINARY_DIR}/conf
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/conf
)

add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/pairinteraction ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/app.py
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/startgui.cmakein ${CMAKE_BINARY_DIR}/pairinteraction
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pairinteraction ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pyqtgraph/pyqtgraph ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/pyqtgraph
  COMMAND ${PYUIC_BINARY} --output ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/plotter.py ${CMAKE_CURRENT_SOURCE_DIR}/plotter.ui
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/adapt.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/pairinteraction ${CMAKE_CURRENT_SOURCE_DIR}/pyqtgraph ${CMAKE_CURRENT_SOURCE_DIR}/plotter.ui ${CMAKE_CURRENT_SOURCE_DIR}/startgui.cmakein
)

add_custom_target(pythonscripts ALL
  DEPENDS ${CMAKE_BINARY_DIR}/pairinteraction ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction
)
add_dependencies(pythonscripts gui)

# Compile Python scripts

if ( CMAKE_HOST_WIN32 OR CMAKE_HOST_APPLE )

  find_program(PYINSTALLER NAMES pyinstaller)

  add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/dist/pairinteraction ${CMAKE_BINARY_DIR}/dist/pairinteraction/conf ${CMAKE_BINARY_DIR}/dist/pairinteraction/databases/quantum_defects.db
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pairinteraction.spec ${CMAKE_BINARY_DIR}/pairinteraction.spec
    COMMAND ${PYINSTALLER} -y ${CMAKE_BINARY_DIR}/pairinteraction.spec
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/conf ${CMAKE_BINARY_DIR}/dist/pairinteraction/conf
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS ${CMAKE_BINARY_DIR}/pairinteraction ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction ${CMAKE_CURRENT_SOURCE_DIR}/conf
  )

  add_custom_target(pyinstaller ALL
    DEPENDS ${CMAKE_BINARY_DIR}/dist/pairinteraction ${CMAKE_BINARY_DIR}/dist/pairinteraction/conf
  )
  add_dependencies(pyinstaller pythonscripts)

  add_dependencies(pyinstaller _pairinteraction_real)
  add_dependencies(pyinstaller _pairinteraction_complex)
  add_dependencies(pyinstaller database)

else( )

  find_package(PythonInterp 3 REQUIRED)
  add_custom_command(TARGET pythonscripts POST_BUILD
    COMMAND ${PYTHON_EXECUTABLE} -m compileall -q ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction
  )

endif( )

# Install Python scripts

if( CMAKE_HOST_APPLE AND WITH_DMG )

  install(DIRECTORY ${CMAKE_BINARY_DIR}/dist/pairinteraction/. DESTINATION gui PATTERN "pairinteraction" EXCLUDE)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/dist/pairinteraction/pairinteraction DESTINATION gui)

  find_program( FILEICON NAMES fileicon )
  if( FILEICON )
    message( STATUS "Found fileicon: ${FILEICON}" )
  else( )
    message( FATAL_ERROR "Could not find fileicon" )
  endif( )

  set(filepath \${CMAKE_INSTALL_PREFIX}/gui/pairinteraction)
  set(iconpath ${CMAKE_MACOSX_GOODIES_PATH}/pairinteraction.icns)
  install(CODE "execute_process(COMMAND ${FILEICON} set ${filepath} ${iconpath})")

elseif( CMAKE_HOST_UNIX )

  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/conf DESTINATION share/pairinteraction/gui)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction DESTINATION share/pairinteraction/gui)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/pairinteraction DESTINATION share/pairinteraction)

  install(CODE "file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tmp)")
  install(CODE "EXECUTE_PROCESS(COMMAND ln -sf ../share/pairinteraction/pairinteraction ${CMAKE_BINARY_DIR}/tmp/pairinteraction)")
  install(PROGRAMS ${CMAKE_BINARY_DIR}/tmp/pairinteraction DESTINATION bin)
  install(CODE "file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/tmp)")

endif( )
