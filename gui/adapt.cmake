file(READ "${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/plotter.py" TMPTXT)
string(REPLACE "from pyqtgraph import" "from .pyqtgraph import" TMPTXT "${TMPTXT}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/plotter.py" "${TMPTXT}")
