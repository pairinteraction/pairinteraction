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

add_custom_target(gui ALL)
add_custom_command(TARGET gui
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/conf ${CMAKE_CURRENT_BINARY_DIR}/conf
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pairinteraction ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pyqtgraph/pyqtgraph ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/pyqtgraph
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/startgui ${CMAKE_CURRENT_BINARY_DIR}/startgui
  COMMAND ${PYUIC} --output ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/plotter.py ${CMAKE_CURRENT_SOURCE_DIR}/plotter.ui
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/adapt.cmake
  COMMAND ${PYTHON} -m compileall ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/app.py
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/startgui ${CMAKE_BINARY_DIR}/pairinteraction
)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/conf DESTINATION share/pairinteraction/gui)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction DESTINATION share/pairinteraction/gui)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/startgui DESTINATION share/pairinteraction/gui)

if(CMAKE_HOST_UNIX)
  install(CODE "file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tmp)")
  install(CODE "EXECUTE_PROCESS(COMMAND ln -sf ../share/pairinteraction/gui/startgui ${CMAKE_BINARY_DIR}/tmp/pairinteraction)")
  install(PROGRAMS ${CMAKE_BINARY_DIR}/tmp/pairinteraction DESTINATION bin)
  install(CODE "file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/tmp)")
endif(CMAKE_HOST_UNIX)
