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

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

find_program( PYUIC NAMES py3uic5 pyuic5 )
if( PYUIC )
  message( STATUS "Found pyuic: ${PYUIC}" )
else( )
  message( FATAL_ERROR "Could not find pyuic" )
endif( )

find_program( PYTHON NAMES python3 )
if( PYTHON )
  message( STATUS "Found python: ${PYTHON}" )
else( )
  message( FATAL_ERROR "Could not find python" )
endif( )

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/conf
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/conf ${CMAKE_CURRENT_BINARY_DIR}/conf
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/conf/*
)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction ${CMAKE_CURRENT_BINARY_DIR}/startgui
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pairinteraction ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pyqtgraph/pyqtgraph ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/pyqtgraph
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/startgui ${CMAKE_CURRENT_BINARY_DIR}/startgui
  COMMAND ${PYUIC} --output ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/plotter.py ${CMAKE_CURRENT_SOURCE_DIR}/plotter.ui
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/adapt.cmake
  COMMAND ${PYTHON} -m compileall ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/pairinteraction/* ${CMAKE_CURRENT_SOURCE_DIR}/pyqtgraph/pyqtgraph/* ${CMAKE_CURRENT_SOURCE_DIR}/startgui ${CMAKE_CURRENT_SOURCE_DIR}/plotter.ui
)

add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/pairinteraction
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/startgui ${CMAKE_BINARY_DIR}/pairinteraction
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/startgui
)

add_custom_target(gui ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction ${CMAKE_CURRENT_BINARY_DIR}/conf ${CMAKE_CURRENT_BINARY_DIR}/startgui ${CMAKE_BINARY_DIR}/pairinteraction
)

if(DMG)

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/dist/startgui
    COMMAND ${CMAKE_COMMAND} -DDMG=${DMG} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/pyinstaller.cmake
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction ${CMAKE_CURRENT_BINARY_DIR}/startgui
  )
  
  add_custom_target(pyinstaller ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/dist/startgui
  )
  
  add_dependencies(pyinstaller gui)

  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/dist/startgui DESTINATION .)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/conf DESTINATION .)

else(DMG)

  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/conf DESTINATION share/pairinteraction/gui)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction DESTINATION share/pairinteraction/gui)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/startgui DESTINATION share/pairinteraction/gui)

  if(CMAKE_HOST_UNIX)
    install(CODE "file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tmp)")
    install(CODE "EXECUTE_PROCESS(COMMAND ln -sf ../share/pairinteraction/gui/startgui ${CMAKE_BINARY_DIR}/tmp/pairinteraction)")
    install(PROGRAMS ${CMAKE_BINARY_DIR}/tmp/pairinteraction DESTINATION bin)
    install(CODE "file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/tmp)")
  endif(CMAKE_HOST_UNIX)

endif(DMG)
