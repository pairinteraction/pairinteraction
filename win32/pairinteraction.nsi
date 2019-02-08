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

    ClearErrors
    ReadRegDword $0 HKLM "SOFTWARE\Microsoft\VisualStudio\14.0\VC\Runtimes\x64" "Installed"
    ReadRegDword $1 HKLM "SOFTWARE\Microsoft\VisualStudio\14.0\VC\Runtimes\x64" "Major"
    ReadRegDword $2 HKLM "SOFTWARE\Microsoft\VisualStudio\14.0\VC\Runtimes\x64" "Minor"
    ReadRegDword $3 HKLM "SOFTWARE\Microsoft\VisualStudio\14.0\VC\Runtimes\x64" "Bld"
    IfErrors update 0

    StrCmp $0 "1" 0 update 
    IntCmp $1 14 0 update 0
    IntCmp $2 0 0 update 0
    IntCmp $3 24215 0 update 0

    Goto next

    update:
      File "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\redist\1033\vcredist_x64.exe"
      ExecWait "$INSTDIR\vcredist_x64.exe"
      Delete "$INSTDIR\vcredist_x64.exe"

      ReadRegDword $0 HKLM "SOFTWARE\Microsoft\VisualStudio\14.0\VC\Runtimes\x64" "Installed"
      StrCmp $0 "1" 0 fail

      Goto next

      fail:
        MessageBox MB_OK "Redistributable for Visual Studio 2015 could not be found!"
    next:
  SectionEnd
SectionGroupEnd

SectionGroup /e "${APP_NAME}"
  Section 'Backend'
    SectionIn RO
    SetOutPath "$INSTDIR"
    File "..\LICENSE*"
    SetOutPath "$INSTDIR\pairinteraction"
    File "${BUILD_DIR}\pairinteraction\Release\*"
    File "${BUILD_DIR}\pairinteraction\pireal.py"
    File "${BUILD_DIR}\pairinteraction\picomplex.py"
    SetOutPath "$INSTDIR\pairinteraction\databases"
    File "${BUILD_DIR}\pairinteraction\databases\*.db"

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
  RMDir /r "$INSTDIR\pairinteraction"
  RMDir /r "$INSTDIR\gui"

  delete "$INSTDIR\LICENSE*"

  delete "$INSTDIR\pairinteraction.ico"
  delete "$INSTDIR\pairinteraction.bat"
  delete "$DESKTOP\pairinteraction.lnk"

  delete "$INSTDIR\uninstall.exe"

  RMDIR "$INSTDIR"

  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APP_NAME}"
SectionEnd
