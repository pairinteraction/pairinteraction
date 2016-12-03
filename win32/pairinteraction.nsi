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

!if ${ARCH} == "x86_64"
  !define ARCHNAME "x86_64"
  !define PROGDIR $PROGRAMFILES64
!else if ${ARCH} == "i686"
  !define ARCHNAME "x86"
  !define PROGDIR $PROGRAMFILES
!else
  !error "Architecture '${ARCH}' unkown"
!endif


!include "MUI.nsh"

!define APP_NAME "pairinteraction"
!define BUILD_DIR "..\build"
!define DLL_DIR "\usr\${ARCH}-w64-mingw32\sys-root\mingw\bin\"
name ${APP_NAME}

OutFile '${APP_NAME}-install-windows-${ARCH}.exe'

showinstdetails show

InstallDir '${PROGDIR}\${APP_NAME}'

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "..\LICENSE"
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_COMPONENTS
!define MUI_FINISHPAGE_NOAUTOCLOSE
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_LANGUAGE "English"

SectionGroup /e "Dependencies"
  Section 'MS-MPI'
    SetOutPath "$INSTDIR"
    File "${BUILD_DIR}\win32\MSMpiSetup.exe"
    ExecWait "$INSTDIR\MSMpiSetup.exe"
    Delete "$INSTDIR\MSMpiSetup.exe"
  SectionEnd

  Section 'Miniconda'
    SetOutPath "$INSTDIR"
    File "${BUILD_DIR}\win32\Miniconda3-latest-Windows-${ARCHNAME}.exe"
    ExecWait "$INSTDIR\Miniconda3-latest-Windows-${ARCHNAME}.exe"

    !define CONDA_PATH "${PROGDIR}\Miniconda3\Scripts\conda.exe"
    !define PIP_PATH "${PROGDIR}\Miniconda3\Scripts\pip.exe"

    IfFileExists "${CONDA_PATH}" 0 fail
    IfFileExists "${PIP_PATH}" 0 fail

    ExecWait "${CONDA_PATH} install -y numpy scipy"
    ExecWait "${PIP_PATH} install psutil pint pyqt5"
    Goto done

fail:
      MessageBox MB_OK "Miniconda installation could not be found!"
done:

    Delete "$INSTDIR\Miniconda3-latest-Windows-${ARCHNAME}.exe"
  SectionEnd
SectionGroupEnd

SectionGroup /e "${APP_NAME}"
  Section 'Backend'
    SectionIn RO
    SetOutPath "$INSTDIR\calc"
    File "${BUILD_DIR}\calc\*.exe"
    File "${DLL_DIR}\libboost_filesystem-mt.dll"
    File "${DLL_DIR}\libboost_system-mt.dll"
    !if ${ARCH} == "x86_64"
      File "${DLL_DIR}\libgcc_s_seh-1.dll"
    !else if ${ARCH} == "i686"
      File "${DLL_DIR}\libgcc_s_sjlj-1.dll"
    !endif
    File "${DLL_DIR}\libgsl-0.dll"
    File "${DLL_DIR}\libgslcblas-0.dll"
    File "${DLL_DIR}\libsqlite3-0.dll"
    File "${DLL_DIR}\libstdc++-6.dll"
    File "${DLL_DIR}\libwinpthread-1.dll"
    SetOutPath "$INSTDIR\calc\databases"
    File "${BUILD_DIR}\calc\databases\*.db"
  SectionEnd

  Section 'GUI (Recommended)'
    SetOutPath "$INSTDIR\gui"
    File "${BUILD_DIR}\gui\startgui"
    File /r "${BUILD_DIR}\gui\conf"
    File /r "${BUILD_DIR}\gui\pairinteraction"
  SectionEnd
SectionGroupEnd

Section 'Desktop Icon'
  !define PYTHON_PATH "C:\Miniconda3\python.exe"
  IfFileExists "${PYTHON_PATH}" success 0
  MessageBox MB_OK "Miniconda installation could not be found!  Assuming python.exe to be in PATH."
  !undef PYTHON_PATH
  !define PYTHON_PATH "python"

success:
  SetOutPath "$INSTDIR\"
  File "pairinteraction.ico"
  FileOpen  $4 "$INSTDIR\pairinteraction.bat" w
  FileWrite $4 "@echo off$\r$\n"
  FileWrite $4 '${PYTHON_PATH} "$INSTDIR\gui\startgui"'
  FileClose $4

  CreateShortCut "$DESKTOP\pairinteraction.lnk" "$INSTDIR\pairinteraction.bat" "" \
    "$INSTDIR\pairinteraction.ico"
SectionEnd
