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

# pip3 install git+https://github.com/pyinstaller/pyinstaller
find_program( PYINSTALLER NAMES pyinstaller )
if( PYINSTALLER )
  message( STATUS "Found pyinstaller: ${PYINSTALLER}" )
else( )
  message( FATAL_ERROR "Could not find pyinstaller" )
endif( )

execute_process(COMMAND "${PYINSTALLER}" "--onefile" "--exclude-module" "matplotlib" "--exclude-module" "OpenGL" "--exclude-module" "PyQt5.QtOpenGL" "${CMAKE_CURRENT_BINARY_DIR}/startgui")
