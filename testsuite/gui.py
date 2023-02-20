# https://johnnado.com/pyqt-qtest-example/
# https://github.com/jmcgeheeiv/pyqttestexample
import sys
import time
import unittest

from PyQt5.QtCore import Qt
from PyQt5.QtTest import QTest
from PyQt5.QtWidgets import QApplication

import pairinteraction_gui.pairinteraction.app as piGui

app = QApplication(sys.argv)


class PairinteractionGuiTest(unittest.TestCase):
    def setUp(self):
        self.form = piGui.MainWindow()
        self.form.ui.action_sconf_reset.trigger()
        self.form.ui.action_pconf_reset.trigger()

    def testFieldmapCalcButton(self):
        self.form.ui.lineedit_system_maxEx.setText("10")
        self.form.ui.spinbox_system_steps.setValue(10)
        widget = self.form.ui.pushbutton_field1_calc
        QTest.mouseClick(widget, Qt.LeftButton)

    def testPotentialCalcButton(self):
        self.form.ui.spinbox_system_steps.setValue(10)
        widget = self.form.ui.pushbutton_potential_calc
        QTest.mouseClick(widget, Qt.LeftButton)

    def tearDown(self):
        # Calculation runs in the background. Wait for it to finish.
        if self.form.thread.isRunning():
            self.form.thread.wait()
        # Close any pipes and wait for subprocess to exit.
        if self.form.proc:
            self.form.proc.stdout.close()
            self.form.proc.wait()


if __name__ == "__main__":
    unittest.main()
