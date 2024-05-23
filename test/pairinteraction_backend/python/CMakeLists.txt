add_test(NAME pytest
  COMMAND ${Python_EXECUTABLE} -m pytest
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(pytest PROPERTIES
  ENVIRONMENT "PYTHONPATH=$<SHELL_PATH:$<TARGET_FILE_DIR:pairinteraction_backend>;$ENV{PYTHONPATH}>")
