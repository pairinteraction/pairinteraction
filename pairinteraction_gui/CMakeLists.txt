# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
#
# This file is part of the pairinteraction GUI.
#
# The pairinteraction GUI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction GUI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the pairinteraction GUI. If not, see <http://www.gnu.org/licenses/>.

find_program(PYUIC_BINARY NAMES py3uic5 pyuic5 py3uic5.bat pyuic5.bat REQUIRED)
if(NOT PYUIC_BINARY)
  message(FATAL_ERROR "Python User Interface Compiler not found")
endif()

set(STARTGUI "start_pairinteraction_gui")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.py.cmakein ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction_app/version.py)
add_custom_target(gui ALL
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/conf ${CMAKE_CURRENT_BINARY_DIR}/conf
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/startgui.cmakein ${CMAKE_BINARY_DIR}/${STARTGUI}
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pairinteraction_app ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction_app
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pipy ${CMAKE_CURRENT_BINARY_DIR}/pipy
  COMMAND ${PYUIC_BINARY} --output ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction_app/plotter.py ${CMAKE_CURRENT_SOURCE_DIR}/plotter.ui
)

# Compile Python scripts

if ( CMAKE_HOST_WIN32 OR (CMAKE_HOST_APPLE AND WITH_DMG) )

  find_program(PYINSTALLER NAMES pyinstaller)

  add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/dist/pairinteraction_gui
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pairinteraction.spec ${CMAKE_BINARY_DIR}/pairinteraction.spec
    COMMAND ${PYINSTALLER} -y ${CMAKE_BINARY_DIR}/pairinteraction.spec
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS ${CMAKE_BINARY_DIR}/${STARTGUI} ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction_app
  )

  add_custom_target(pyinstaller ALL
    DEPENDS ${CMAKE_BINARY_DIR}/dist/pairinteraction_gui
  )
  add_dependencies(pyinstaller gui)

  add_dependencies(pyinstaller _binding)
  add_dependencies(pyinstaller database)

else( )

  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  add_custom_command(TARGET gui POST_BUILD
    COMMAND ${Python3_EXECUTABLE} -m compileall -q ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction_app
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction_app
  )

endif( )

# Install Python scripts

if( CMAKE_HOST_APPLE AND WITH_DMG )

  install(DIRECTORY ${CMAKE_BINARY_DIR}/dist/pairinteraction_gui/. DESTINATION pairinteraction_gui PATTERN "pairinteraction_gui" EXCLUDE)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/dist/pairinteraction_gui/pairinteraction_gui DESTINATION pairinteraction_gui)

  find_program( FILEICON NAMES fileicon )
  if( FILEICON )
    message( STATUS "Found fileicon: ${FILEICON}" )
  else( )
    message( FATAL_ERROR "Could not find fileicon" )
  endif( )

  set(filepath \${CMAKE_INSTALL_PREFIX}/pairinteraction_gui/pairinteraction_gui)
  set(iconpath ${CMAKE_MACOSX_GOODIES_PATH}/pairinteraction.icns)
  install(CODE "execute_process(COMMAND ${FILEICON} set ${filepath} ${iconpath})")

elseif( CMAKE_HOST_UNIX )

  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/conf DESTINATION share/pairinteraction/pairinteraction_gui)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction_app DESTINATION share/pairinteraction/pairinteraction_gui)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pipy DESTINATION share/pairinteraction/pairinteraction_gui)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/${STARTGUI} DESTINATION share/pairinteraction)

  install(CODE "file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tmp)")
  install(CODE "EXECUTE_PROCESS(COMMAND ln -sf ../share/pairinteraction/${STARTGUI} ${CMAKE_BINARY_DIR}/tmp/${STARTGUI})")
  install(PROGRAMS ${CMAKE_BINARY_DIR}/tmp/${STARTGUI} DESTINATION bin)
  install(CODE "file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/tmp)")

endif( )
