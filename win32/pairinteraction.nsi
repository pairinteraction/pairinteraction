!include "MUI.nsh"

!define APP_NAME "pairinteraction"
!define BUILD_DIR "..\build-w64"
name ${APP_NAME}

OutFile '${APP_NAME}-install-x86_64.exe'

showinstdetails show

InstallDir '$PROGRAMFILES64\${APP_NAME}'

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "license.txt"
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_COMPONENTS
!define MUI_FINISHPAGE_NOAUTOCLOSE
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_LANGUAGE "English"

SectionGroup /e "Dependencies"
  Section 'MS-MPI'
    SetOutPath "$INSTDIR"
    File "MSMpiSetup.exe"
    ExecWait "$INSTDIR\MSMpiSetup.exe"
    Delete "$INSTDIR\MSMpiSetup.exe"
  SectionEnd

  Section 'Miniconda'
    SetOutPath "$INSTDIR"
    File "Miniconda3-latest-Windows-x86_64.exe"
    ExecWait "$INSTDIR\Miniconda3-latest-Windows-x86_64.exe"

    !define CONDA_PATH "C:\Miniconda3\Scripts\conda.exe"
    !define PIP_PATH "C:\Miniconda3\Scripts\pip.exe"

    IfFileExists "${CONDA_PATH}" 0 fail
    IfFileExists "${PIP_PATH}" 0 fail

    ExecWait "${CONDA_PATH} install -y numpy scipy"
    ExecWait "${PIP_PATH} install psutil pint pyqt5"
    Goto done

    fail:
      MessageBox MB_OK "Miniconda installation could not be found!"
    done:

    Delete "$INSTDIR\Miniconda3-latest-Windows-x86_64.exe"
  SectionEnd
SectionGroupEnd

SectionGroup /e "${APP_NAME}"
  Section 'Backend'
    SetOutPath "$INSTDIR\calc"
    File "${BUILD_DIR}\calc\*.exe"
    File "${BUILD_DIR}\calc\*.dll"
    File "${BUILD_DIR}\calc\wignerSymbols\*.dll"
    File "${BUILD_DIR}\calc\databases\*.db"
  SectionEnd

  Section 'GUI'
    SetOutPath "$INSTDIR\gui"
    File "${BUILD_DIR}\gui\*.py"
    File "${BUILD_DIR}\gui\*conf"
    File /r "${BUILD_DIR}\gui\__pycache__"
    File /r "${BUILD_DIR}\gui\pyqtgraph"
  SectionEnd
SectionGroupEnd

Section 'Desktop Icon'
  !define PYTHON_PATH "C:\Miniconda3\python.exe"
  IfFileExists "${PYTHON_PATH}" +4 0
  MessageBox MB_OK "Miniconda installation could not be found!  Assuming python.exe to be in PATH."
  !undef PYTHON_PATH
  !define PYTHON_PATH "python"

  FileOpen  $4 "$INSTDIR\pairinteraction.bat" w
  FileWrite $4 "@echo off"
  FileWrite $4 '${PYTHON_PATH} "$INSTDIR\gui\pairinteraction.py"'
  FileClose $4

  CreateShortCut "$DESKTOP\pairinteraction.lnk" "$INSTDIR\pairinteraction.bat"
SectionEnd
