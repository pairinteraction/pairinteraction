cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

add_custom_target(gui ALL)
add_custom_command(TARGET gui
  COMMAND pyuic4 ${CMAKE_CURRENT_SOURCE_DIR}/plotter.ui > ${CMAKE_CURRENT_BINARY_DIR}/plotter.py
  COMMAND python -m compileall ${CMAKE_CURRENT_BINARY_DIR}/plotter.py
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/start.py ${CMAKE_CURRENT_BINARY_DIR}/start.py
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/start.py ${CMAKE_BINARY_DIR}/pairinteraction-gui
  COMMAND ${CMAKE_COMMAND} -E create_symlink ../share/pairinteraction/gui/start.py ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction-gui
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/plotter.py DESTINATION share/pairinteraction/gui)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/plotter.pyc DESTINATION share/pairinteraction/gui)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/start.py DESTINATION share/pairinteraction/gui)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pairinteraction-gui DESTINATION bin)
