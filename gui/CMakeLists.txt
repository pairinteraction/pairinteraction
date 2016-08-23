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
  COMMAND ${PYUIC} ${CMAKE_CURRENT_SOURCE_DIR}/plotter.ui > ${CMAKE_CURRENT_BINARY_DIR}/plotter.py
  COMMAND ${PYTHON} -m compileall ${CMAKE_CURRENT_BINARY_DIR}/plotter.py
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/pairinteraction.py ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction.py
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/example.sconf ${CMAKE_CURRENT_BINARY_DIR}/example.sconf
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/example.pconf ${CMAKE_CURRENT_BINARY_DIR}/example.pconf
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pyqtgraph/pyqtgraph ${CMAKE_CURRENT_BINARY_DIR}/pyqtgraph
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction.py ${CMAKE_BINARY_DIR}/pairinteraction
  COMMAND ${CMAKE_COMMAND} -E create_symlink ../share/pairinteraction/gui/pairinteraction.py ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/example.sconf DESTINATION share/pairinteraction/gui)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/example.pconf DESTINATION share/pairinteraction/gui)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/plotter.py DESTINATION share/pairinteraction/gui)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/__pycache__ DESTINATION share/pairinteraction/gui)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pyqtgraph DESTINATION share/pairinteraction/gui)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction.py DESTINATION share/pairinteraction/gui)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction DESTINATION bin)
