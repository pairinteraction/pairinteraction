#!/usr/bin/env bash

# http://stackoverflow.com/questions/26465240/how-to-execute-a-shell-script-from-within-an-app-and-show-it-in-the-terminal

scriptPath=$(dirname "$0")'/../Resources/pairinteraction/gui/start_pairinteraction_gui.py'
osascript -e 'tell app "Terminal" to do script "\"'"$scriptPath"'\"; exit"'
