install(FILES metainfo.xml DESTINATION share/metainfo/$ENV{FLATPAK_ID}.metainfo.xml)
install(FILES pairinteraction.desktop DESTINATION share/applications/$ENV{FLATPAK_ID}.desktop)
install(FILES pairinteraction.png DESTINATION share/icons/hicolor/512x512/apps/$ENV{FLATPAK_ID}.png)
