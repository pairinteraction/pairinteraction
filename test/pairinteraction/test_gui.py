# https://johnnado.com/pyqt-qtest-example/
# https://github.com/jmcgeheeiv/pyqttestexample
import io
import json
import os
import shutil
import sys
import tempfile
import unittest
import zipfile

import numpy as np
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtTest import QTest
from scipy.io import loadmat

from pairinteraction_gui.app import MainWindow

PATH = os.path.join("data")

app = QtWidgets.QApplication(sys.argv)


class PairinteractionGuiTest(unittest.TestCase):
    def setUp(self):
        self.form = MainWindow()
        self.form.path_cache = tempfile.mkdtemp()
        self.form.ui.action_sconf_reset.trigger()
        self.form.ui.action_pconf_reset.trigger()

    def testFieldCalcButton(self):
        # Testing simulation single atom with all E and B fields on
        self.form.loadSettingsSystem(os.path.join(PATH, "field", "settings.sconf"))
        self.form.loadSettingsPlotter(os.path.join(PATH, "field", "settings.pconf"))
        self._testEnergies(0, "field", dE=3)

    def testPotentialCalcButton(self):
        # Testing simulation for pairpotential
        self.form.loadSettingsSystem(os.path.join(PATH, "potential", "settings.sconf"))
        self.form.loadSettingsPlotter(os.path.join(PATH, "potential", "settings.pconf"))
        self._testEnergies(2, "potential", dE=0.3)

    @unittest.skip("TODO implement new simulation call")
    def _testEnergies(self, idx, ref_data, dE, dE_tol=1e-3, use_python_api="default"):
        if use_python_api == "default":
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
        self.form.savePlot = False
        forceFilename = os.path.join(PATH, "current")
        self.form.forceFilename = forceFilename
        QTest.mouseClick(widget_save, Qt.LeftButton)

        data = {}
        sconfig = {}
        # Load current data
        with zipfile.ZipFile(forceFilename, "r") as zip_file:
            with zip_file.open("data.mat") as f:
                f_io = io.BytesIO(f.read())
                data["current"] = loadmat(f_io)
            with zip_file.open("settings.sconf") as f:
                sconfig["current"] = json.load(f)
        os.remove(forceFilename)

        # Load reference data
        with open(os.path.join(PATH, ref_data, "data.mat"), "rb") as f:
            data["ref"] = loadmat(f)
        with open(os.path.join(PATH, ref_data, "settings.sconf")) as f:
            sconfig["ref"] = json.load(f)

        # Check if configs match # unecessary since we load the same config
        for k, v in sconfig["ref"].items():
            assert sconfig["current"][k] == v

        if len(data["current"]["eigenvalues"]) == 0:
            raise ValueError("No eigenvalues found in current data, some bug in calculating the eigenvalues.")

        # Check if central eigenvalues (+/- dE) match
        for i in range(len(data["ref"]["eigenvalues"])):
            Es = {k: np.array(mat["eigenvalues"])[i] for k, mat in data.items()}
            Es = {k: E[np.abs(E) < dE] for k, E in Es.items()}
            assert len(np.abs(Es["ref"]) == len(Es["current"]))
            diff_rel = np.abs(Es["ref"] - Es["current"])
            assert np.all(diff_rel <= dE_tol)

    def tearDown(self):
        # clean up processes
        self.form.cleanupProcesses()

        # Remove tmp cache
        shutil.rmtree(self.form.path_cache)


if __name__ == "__main__":
    unittest.main()
