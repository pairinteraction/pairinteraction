cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

FIND_PROGRAM(PYUIC NAMES py3uic5 pyuic5)
IF(PYUIC)
  MESSAGE(STATUS "Found pyuic: ${PYUIC}")
ELSE(PYUIC)
  MESSAGE(FATAL_ERROR "Could not find pyuic")
ENDIF(PYUIC)

FIND_PROGRAM(PYTHON NAMES python3)
IF(PYUIC)
  MESSAGE(STATUS "Found python: ${PYTHON}")
ELSE(PYUIC)
  MESSAGE(FATAL_ERROR "Could not find python")
ENDIF(PYUIC)

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
