# https://johnnado.com/pyqt-qtest-example/
# https://github.com/jmcgeheeiv/pyqttestexample
import json
import os
import sys
import unittest
import zipfile

import numpy as np
import scipy.io
from PyQt5.QtCore import Qt
from PyQt5.QtTest import QTest
from PyQt5.QtWidgets import QApplication

import pairinteraction_gui.pairinteraction.app as piGui

app = QApplication(sys.argv)
PATH = "reference_data/gui/"


class PairinteractionGuiTest(unittest.TestCase):
    def setUp(self):
        self.form = piGui.MainWindow()
        self.form.ui.action_sconf_reset.trigger()
        self.form.ui.action_pconf_reset.trigger()

    def testFieldCalcButton(self):
        # Testing simulation single atom with all E and B fields on
        self.form.loadSettingsSystem(PATH + "Field/settings.sconf")
        self.form.loadSettingsPlotter(PATH + "Field/settings.pconf")
        self._testEnergies(0, "Field", dE=3)

    def testPotentialCalcButton(self):
        # Testing simulation for pairpotential
        self.form.loadSettingsSystem(PATH + "Potential/settings.sconf")
        self.form.loadSettingsPlotter(PATH + "Potential/settings.pconf")
        self._testEnergies(2, "Potential", dE=0.3)

    def _testEnergies(self, idx, ref_data, dE, dE_tol=1e-3, use_python_api="both"):
        if use_python_api == "both":
            for use_python_api in [False, True]:
                self._testEnergies(idx, ref_data, dE, use_python_api=use_python_api)
            return
        self.form.ui.checkbox_use_python_api.setChecked(use_python_api)
        self.form.autosetSymmetrization()

        if idx == 0:
            widget_calc = self.form.ui.pushbutton_field1_calc
            widget_save = self.form.ui.pushbutton_field1_save
        elif idx == 2:
            widget_calc = self.form.ui.pushbutton_potential_calc
            widget_save = self.form.ui.pushbutton_potential_save

        # Run simulation
        QTest.mouseClick(widget_calc, Qt.LeftButton)
        while self.form.timer.isActive():
            self.form.checkForData()

        # Save current data
        self.form.forceFilename = PATH + "tmp"
        QTest.mouseClick(widget_save, Qt.LeftButton)

        data = {}
        sconfig = {}
        # Load reference data
        with zipfile.ZipFile(PATH + "tmp", "r") as zip_file:
            with zip_file.open("data.mat") as f:
                data["tmp"] = scipy.io.loadmat(f)
            with zip_file.open("settings.sconf") as f:
                sconfig["tmp"] = json.load(f)
        os.remove(PATH + "tmp")

        # Load current data
        with open(PATH + ref_data + "/data.mat", "rb") as f:
            data["ref"] = scipy.io.loadmat(f)
        with open(PATH + ref_data + "/settings.sconf") as f:
            sconfig["ref"] = json.load(f)

        # Check if configs match # unecessary since we load the same config
        for k, v in sconfig["ref"].items():
            assert sconfig["tmp"][k] == v

        # Check if central eigenvalues (+/- dE) match
        for i in range(len(data["ref"]["eigenvalues"])):
            Es = {k: np.array(mat["eigenvalues"])[i] for k, mat in data.items()}
            Es = {k: E[np.abs(E) < dE] for k, E in Es.items()}
            diff_rel = np.abs(Es["ref"] - Es["tmp"])
            assert np.all(diff_rel <= dE_tol)

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
