# https://johnnado.com/pyqt-qtest-example/
# https://github.com/jmcgeheeiv/pyqttestexample

import time
import sys
import unittest
from PyQt5.QtWidgets import QApplication
from PyQt5.QtTest import QTest
from PyQt5.QtCore import Qt

import pairinteraction_gui.pairinteraction.app as piGui

app = QApplication(sys.argv)


class PairinteractionGuiTest(unittest.TestCase):

    def setUp(self):
        self.form = piGui.MainWindow()

    def testPotentialCalcButton(self):
        widget = self.form.ui.pushbutton_potential_calc
        QTest.mouseClick(widget, Qt.LeftButton)

        # Calculation runs in the background. Wait for it to finish.
        self.form.thread.wait()

        # Close any pipes and wait for subprocess to exit.
        self.form.proc.stdout.close()
        self.form.proc.wait()


if __name__ == "__main__":
    unittest.main()
