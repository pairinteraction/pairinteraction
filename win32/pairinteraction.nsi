# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
#
# This file is part of the pairinteraction library.
#
# The pairinteraction library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.

!include "MUI.nsh"

!define APP_NAME "pairinteraction"
!define LIBNAME "pairinteraction"
!define GUINAME "pairinteraction_gui"
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

    !getdllversion "${BUILD_DIR}\vcredist_x64.exe" RedistVer

    StrCmp $0 "1" 0 update
    IntCmp $1 "${RedistVer1}" 0 update 0
    IntCmp $2 "${RedistVer2}" 0 update 0
    IntCmp $3 "${RedistVer3}" 0 update 0

    Goto next

    update:
      File "${BUILD_DIR}\vcredist_x64.exe"
      ExecWait "$INSTDIR\vcredist_x64.exe"
      Delete "$INSTDIR\vcredist_x64.exe"

      ReadRegDword $0 HKLM "SOFTWARE\Microsoft\VisualStudio\14.0\VC\Runtimes\x64" "Installed"
      StrCmp $0 "1" 0 fail

      Goto next

      fail:
        MessageBox MB_OK "Redistributable for Visual Studio 2017 could not be found!"
    next:
  SectionEnd
SectionGroupEnd

SectionGroup /e "${APP_NAME}"
  Section 'Backend'
    SectionIn RO
    SetOutPath "$INSTDIR"
    File "..\LICENSE*"
    SetOutPath "$INSTDIR\${LIBNAME}"
    File "${BUILD_DIR}\${LIBNAME}\Release\*"
    File "${BUILD_DIR}\${LIBNAME}\pireal.py"
    File "${BUILD_DIR}\${LIBNAME}\picomplex.py"
    SetOutPath "$INSTDIR\${LIBNAME}\databases"
    File "${BUILD_DIR}\${LIBNAME}\databases\*.db"

    writeUninstaller "$INSTDIR\uninstall.exe"

    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APP_NAME}" \
                     "DisplayName" "${APP_NAME}"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APP_NAME}" \
                     "UninstallString" "$\"$INSTDIR\uninstall.exe$\""
  SectionEnd

  Section 'GUI (Recommended)'
    SetOutPath "$INSTDIR\${GUINAME}"
    File /r "${BUILD_DIR}\dist\pairinteraction\*"
  SectionEnd
SectionGroupEnd


Section 'Desktop Icon'
  SetOutPath "$INSTDIR\"
  File "pairinteraction.ico"
  FileOpen  $4 "$INSTDIR\pairinteraction.bat" w
  FileWrite $4 "@echo off$\r$\n"
  FileWrite $4 'cmd /k ""$INSTDIR\${GUINAME}\pairinteraction.exe""'
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
  RMDir /r "$INSTDIR\${LIBNAME}"
  RMDir /r "$INSTDIR\${GUINAME}"

  delete "$INSTDIR\LICENSE*"

  delete "$INSTDIR\pairinteraction.ico"
  delete "$INSTDIR\pairinteraction.bat"
  delete "$DESKTOP\pairinteraction.lnk"

  delete "$INSTDIR\uninstall.exe"

  RMDIR "$INSTDIR"

  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APP_NAME}"
SectionEnd
