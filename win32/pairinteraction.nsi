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

!include "MUI.nsh"

!define APP_NAME "pairinteraction"
!define BUILD_DIR "..\build"
!define DLL_DIR "..\vcpkg-export\installed\x64-windows\bin\"
!define PROGDIR $PROGRAMFILES64

name ${APP_NAME}

OutFile '${BUILD_DIR}\${APP_NAME}-install-windows.exe'

showinstdetails show

InstallDir '${PROGDIR}\${APP_NAME}'

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "${BUILD_DIR}\LICENSES.txt"
!insertmacro MUI_PAGE_LICENSE "${BUILD_DIR}\..\LICENSE-MKL.txt"
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_COMPONENTS
!define MUI_FINISHPAGE_NOAUTOCLOSE
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_LANGUAGE "English"

SectionGroup /e "Dependencies"
  Section 'Redistributable for Visual Studio 2015'
    SetOutPath "$INSTDIR"
    File "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\redist\1033\vcredist_x64.exe"
    Exec "$INSTDIR\vcredist_x64.exe"

    #ReadRegDword $0 HKLM "SOFTWARE\Wow6432Node\Microsoft\VisualStudio\14.0\VC\Runtimes\x86" "Installed"
    #${If} $0 == "1"
    #      Goto done
    #${Else}
    #      Goto fail
    #${EndIf} # TODO make this check work, see https://stackoverflow.com/questions/12206314/detect-if-visual-c-redistributable-for-visual-studio-2012-is-installed
    Goto done
fail:
      MessageBox MB_OK "Redistributable for Visual Studio 2015 could not be found!"
done:
    Delete "$INSTDIR\vcredist_x64.exe"
  SectionEnd
SectionGroupEnd

SectionGroup /e "${APP_NAME}"
  Section 'Backend'
    SectionIn RO
    SetOutPath "$INSTDIR"
    File "..\LICENSE*"
    SetOutPath "$INSTDIR\libpairinteraction"
    File "${BUILD_DIR}\libpairinteraction\Release\*"
    File "${BUILD_DIR}\libpairinteraction\pireal.py"
    File "${BUILD_DIR}\libpairinteraction\picomplex.py"
    File "${DLL_DIR}\gsl.dll"
    File "${DLL_DIR}\gslcblas.dll"
    File "${DLL_DIR}\libzmq.dll"
    File "${DLL_DIR}\sqlite3.dll"
    File "${DLL_DIR}\boost_filesystem-vc140-mt-1_64.dll"
    File "${DLL_DIR}\boost_program_options-vc140-mt-1_64.dll"
    File "${DLL_DIR}\boost_serialization-vc140-mt-1_64.dll"
    File "${DLL_DIR}\boost_system-vc140-mt-1_64.dll"
    SetOutPath "$INSTDIR\libpairinteraction\databases"
    File "${BUILD_DIR}\libpairinteraction\databases\*.db"

    writeUninstaller "$INSTDIR\uninstall.exe"

    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APP_NAME}" \
                     "DisplayName" "${APP_NAME}"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APP_NAME}" \
                     "UninstallString" "$\"$INSTDIR\uninstall.exe$\""
  SectionEnd

  Section 'GUI (Recommended)'
    SetOutPath "$INSTDIR\gui"
    File /r "${BUILD_DIR}\dist\pairinteraction\*"
  SectionEnd
SectionGroupEnd


Section 'Desktop Icon'
  SetOutPath "$INSTDIR\"
  File "pairinteraction.ico"
  FileOpen  $4 "$INSTDIR\pairinteraction.bat" w
  FileWrite $4 "@echo off$\r$\n"
  FileWrite $4 'cmd /k ""$INSTDIR\gui\pairinteraction.exe""'
  FileClose $4

  CreateShortCut "$DESKTOP\pairinteraction.lnk" "$INSTDIR\pairinteraction.bat" "" \
    "$INSTDIR\pairinteraction.ico"
SectionEnd


function un.onInit
  SetShellVarContext all

  MessageBox MB_OKCANCEL "Do you want to remove ${APP_NAME} and all of its components?" IDOK next
    Abort
  next:
functionEnd


Section 'uninstall'
  RMDir /r "$INSTDIR\libpairinteraction"
  RMDir /r "$INSTDIR\gui"

  delete "$INSTDIR\LICENSE*"

  delete "$INSTDIR\pairinteraction.ico"
  delete "$INSTDIR\pairinteraction.bat"
  delete "$DESKTOP\pairinteraction.lnk"

  delete "$INSTDIR\uninstall.exe"

  RMDIR "$INSTDIR"

  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APP_NAME}"
SectionEnd