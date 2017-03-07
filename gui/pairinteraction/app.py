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

# Standard library
from datetime import timedelta, datetime
from io import StringIO, BytesIO
import locale
import multiprocessing
from operator import itemgetter
import json
import os
import shutil
import signal
import sys
import multiprocessing
from time import sleep, time, strftime
import zipfile

# Communication
import zmq

# Process information
import psutil

# Units
from pint import UnitRegistry
from pint.unit import UndefinedUnitError

# GUI
import sip
from PyQt5 import QtCore, QtGui
from PyQt5.QtPrintSupport import QPrintDialog, QPrinter
from . import pyqtgraph as pg
from .pyqtgraph import exporters
from .pyqtgraph import ColorButton, GradientWidget, PlotWidget
from .plotter import Ui_plotwindow

# Numerics
import numpy as np
from scipy import sparse
from scipy.ndimage.filters import gaussian_filter
from scipy import io

# Own classes
from .utils import Wignerd, csc_happend, csr_vappend, csr_keepmax, bytescale
from .unitmanagement import Quantity, Units
from .guiadditions import GuiDict, DoubledeltaValidator, DoublenoneValidator, DoublepositiveValidator, DoubleValidator
from .pyqtgraphadditions import PointsItem, MultiLine
from .worker import Worker
from .loader import Eigensystem

from calc import pairinteraction_real as pir
from calc import pairinteraction_complex as pic
pi = None

# Versioning
version_settings = 14
version_cache = 12


# Make program killable via strg-c if it is started in a terminal
signal.signal(signal.SIGINT, signal.SIG_DFL)

# Global configurations of pyqtgraph
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')


# === Dictionary to manage the elements of the GUI related to the plotter ===

class PlotDict(GuiDict):

    def _setup(self, store, ui):
        store["minE_field1"] = {
            'widget': ui.lineedit_field1_minE,
            'unit': Units.energy}
        store["maxE_field1"] = {
            'widget': ui.lineedit_field1_maxE,
            'unit': Units.energy}
        store["minE_field2"] = {
            'widget': ui.lineedit_field2_minE,
            'unit': Units.energy}
        store["maxE_field2"] = {
            'widget': ui.lineedit_field2_maxE,
            'unit': Units.energy}
        store["minE_potential"] = {
            'widget': ui.lineedit_potential_minE,
            'unit': Units.energy}
        store["maxE_potential"] = {
            'widget': ui.lineedit_potential_maxE,
            'unit': Units.energy}
        store["lines"] = {
            'widget': ui.groupbox_plot_lines,
            'unit': None}
        store["points"] = {
            'widget': ui.groupbox_plot_points,
            'unit': None}
        store["labels"] = {
            'widget': ui.groupbox_plot_labels,
            'unit': None}
        store["overlap"] = {
            'widget': ui.groupbox_plot_overlap,
            'unit': None}
        store["szLine"] = {
            'widget': ui.spinbox_plot_szLine,
            'unit': None}
        store["szPoint"] = {
            'widget': ui.spinbox_plot_szPoint,
            'unit': None}
        store["szLabel"] = {
            'widget': ui.spinbox_plot_szLabel,
            'unit': None}
        store["szOverlap"] = {
            'widget': ui.spinbox_plot_szOverlap,
            'unit': None}
        store["transpLine"] = {
            'widget': ui.spinbox_plot_transpLine,
            'unit': None}
        store["transpPoint"] = {
            'widget': ui.spinbox_plot_transpPoint,
            'unit': None}
        store["transpLabel"] = {
            'widget': ui.spinbox_plot_transpLabel,
            'unit': None}
        store["transpOverlap"] = {
            'widget': ui.spinbox_plot_transpOverlap,
            'unit': None}
        store["overlapUnperturbed"] = {
            'widget': ui.radiobutton_plot_overlapUnperturbed,
            'unit': None}
        store["overlapDefined"] = {
            'widget': ui.radiobutton_plot_overlapDefined,
            'unit': None}
        store["connectionthreshold"] = {
            'widget': ui.spinbox_plot_connectionthreshold,
            'unit': None}
        store["lin"] = {
            'widget': ui.radiobutton_plot_lin,
            'unit': None}
        store["log"] = {
            'widget': ui.radiobutton_plot_log,
            'unit': None}
        store["resolution"] = {
            'widget': ui.spinbox_plot_resolution,
            'unit': None}
        store["n1"] = {
            'widget': ui.spinbox_plot_n1,
            'unit': None}
        store["n2"] = {
            'widget': ui.spinbox_plot_n2,
            'unit': None}
        store["l1"] = {
            'widget': ui.spinbox_plot_l1,
            'unit': None}
        store["l2"] = {
            'widget': ui.spinbox_plot_l2,
            'unit': None}
        store["j1"] = {
            'widget': ui.spinbox_plot_j1,
            'unit': None}
        store["j2"] = {
            'widget': ui.spinbox_plot_j2,
            'unit': None}
        store["m1"] = {
            'widget': ui.spinbox_plot_m1,
            'unit': None}
        store["m2"] = {
            'widget': ui.spinbox_plot_m2,
            'unit': None}
        store["antialiasing"] = {
            'widget': ui.checkbox_plot_antialiasing,
            'unit': None}
        store["autorange"] = {
            'widget': ui.checkbox_plot_autorange,
            'unit': None}


# === dictionary to manage the elements of the gui related to the system ===

class SystemDict(GuiDict):

    def _setup(self, store, ui):
        store["species1"] = {
            'widget': ui.combobox_system_species1,
            'unit': None}
        store["species2"] = {
            'widget': ui.combobox_system_species2,
            'unit': None}
        store["n1"] = {
            'widget': ui.spinbox_system_n1,
            'unit': None}
        store["n2"] = {
            'widget': ui.spinbox_system_n2,
            'unit': None}
        store["l1"] = {
            'widget': ui.spinbox_system_l1,
            'unit': None}
        store["l2"] = {
            'widget': ui.spinbox_system_l2,
            'unit': None}
        store["j1"] = {
            'widget': ui.spinbox_system_j1,
            'unit': None}
        store["j2"] = {
            'widget': ui.spinbox_system_j2,
            'unit': None}
        store["m1"] = {
            'widget': ui.spinbox_system_m1,
            'unit': None}
        store["m2"] = {
            'widget': ui.spinbox_system_m2,
            'unit': None}
        store["deltaNSingle"] = {
            'widget': ui.spinbox_system_deltaNSingle,
            'unit': None}
        store["deltaLSingle"] = {
            'widget': ui.spinbox_system_deltaLSingle,
            'unit': None}
        store["deltaJSingle"] = {
            'widget': ui.spinbox_system_deltaJSingle,
            'unit': None}
        store["deltaMSingle"] = {
            'widget': ui.spinbox_system_deltaMSingle,
            'unit': None}
        store["deltaNPair"] = {
            'widget': ui.spinbox_system_deltaNPair,
            'unit': None}
        store["deltaLPair"] = {
            'widget': ui.spinbox_system_deltaLPair,
            'unit': None}
        store["deltaJPair"] = {
            'widget': ui.spinbox_system_deltaJPair,
            'unit': None}
        store["deltaMPair"] = {
            'widget': ui.spinbox_system_deltaMPair,
            'unit': None}
        store["deltaESingle"] = {
            'widget': ui.lineedit_system_deltaESingle,
            'unit': Units.energy}
        store["deltaEPair"] = {
            'widget': ui.lineedit_system_deltaEPair,
            'unit': Units.energy}
        store["pairbasisSame"] = {
            'widget': ui.radiobutton_system_pairbasisSame,
            'unit': None}
        store["pairbasisDefined"] = {
            'widget': ui.radiobutton_system_pairbasisDefined,
            'unit': None}
        store["quantizationZ"] = {
            'widget': ui.radiobutton_system_quantizationZ,
            'unit': None}
        store["quantizationInteratomic"] = {
            'widget': ui.radiobutton_system_quantizationInteratomic,
            'unit': None}
        store["samebasis"] = {
            'widget': ui.checkbox_system_samebasis,
            'unit': None}
        store["minEx"] = {
            'widget': ui.lineedit_system_minEx,
            'unit': Units.efield}
        store["minEy"] = {
            'widget': ui.lineedit_system_minEy,
            'unit': Units.efield}
        store["minEz"] = {
            'widget': ui.lineedit_system_minEz,
            'unit': Units.efield}
        store["minBx"] = {
            'widget': ui.lineedit_system_minBx,
            'unit': Units.bfield}
        store["minBy"] = {
            'widget': ui.lineedit_system_minBy,
            'unit': Units.bfield}
        store["minBz"] = {
            'widget': ui.lineedit_system_minBz,
            'unit': Units.bfield}
        store["maxEx"] = {
            'widget': ui.lineedit_system_maxEx,
            'unit': Units.efield}
        store["maxEy"] = {
            'widget': ui.lineedit_system_maxEy,
            'unit': Units.efield}
        store["maxEz"] = {
            'widget': ui.lineedit_system_maxEz,
            'unit': Units.efield}
        store["maxBx"] = {
            'widget': ui.lineedit_system_maxBx,
            'unit': Units.bfield}
        store["maxBy"] = {
            'widget': ui.lineedit_system_maxBy,
            'unit': Units.bfield}
        store["maxBz"] = {
            'widget': ui.lineedit_system_maxBz,
            'unit': Units.bfield}
        store["minR"] = {
            'widget': ui.lineedit_system_minR,
            'unit': Units.length}
        store["maxR"] = {
            'widget': ui.lineedit_system_maxR,
            'unit': Units.length}
        store["theta"] = {
            'widget': ui.lineedit_system_theta,
            'unit': Units.angle}
        store["exponent"] = {
            'widget': ui.spinbox_system_exponent,
            'unit': None}
        store["steps"] = {
            'widget': ui.spinbox_system_steps,
            'unit': None}
        store["precision"] = {
            'widget': ui.lineedit_system_precision,
            'unit': None}
        store["matCombined"] = {
            'widget': ui.radiobutton_system_matCombined,
            'unit': None}
        store["matSeparate"] = {
            'widget': ui.radiobutton_system_matSeparate,
            'unit': None}
        store["missingCalc"] = {
            'widget': ui.radiobutton_system_missingCalc,
            'unit': None}
        store["missingWhittaker"] = {
            'widget': ui.radiobutton_system_missingWhittaker,
            'unit': None}
        store["cores"] = {
            'widget': ui.spinbox_system_cores,
            'unit': None}
        store["diamagnetism"] = {
            'widget': ui.checkbox_system_diamagnetic,
            'unit': None}
        store["symAuto"] = {
            'widget': ui.radiobutton_symAuto,
            'unit': None}
        store["symManual"] = {
            'widget': ui.radiobutton_symManual,
            'unit': None}
        store["invE"] = {
            'widget': ui.checkbox_system_invE,
            'unit': None}
        store["invO"] = {
            'widget': ui.checkbox_system_invO,
            'unit': None}
        store["perE"] = {
            'widget': ui.checkbox_system_perE,
            'unit': None}
        store["perO"] = {
            'widget': ui.checkbox_system_perO,
            'unit': None}
        store["refE"] = {
            'widget': ui.checkbox_system_refE,
            'unit': None}
        store["refO"] = {
            'widget': ui.checkbox_system_refO,
            'unit': None}
        store["conserveM"] = {
            'widget': ui.checkbox_system_conserveM,
            'unit': None}
        store["sametrafo"] = {
            'widget': ui.checkbox_system_sametrafo,
            'unit': None}

    # field map of atom 1 (samebasis == False)
    keys_for_cprogram_field1 = ["species1", "n1", "l1", "j1", "m1",
                                "deltaESingle", "deltaLSingle", "deltaJSingle", "deltaMSingle", "deltaNSingle",
                                "samebasis", "steps", "precision", "missingCalc", "missingWhittaker",
                                "minEx", "minEy", "minEz", "minBx", "minBy", "minBz", "maxEx", "maxEy", "maxEz", "maxBx", "maxBy", "maxBz", "diamagnetism"]

    # field map of atom 2 (samebasis == False)
    keys_for_cprogram_field2 = ["species2", "n2", "l2", "j2", "m2",
                                "deltaESingle", "deltaLSingle", "deltaJSingle", "deltaMSingle", "deltaNSingle",
                                "samebasis", "steps", "precision", "missingCalc", "missingWhittaker",
                                "minEx", "minEy", "minEz", "minBx", "minBy", "minBz", "maxEx", "maxEy", "maxEz", "maxBx", "maxBy", "maxBz", "diamagnetism"]

    # pair potential
    keys_for_cprogram_potential = ["species1", "n1", "l1", "j1", "m1", "species2", "n2", "l2", "j2", "m2",
                                   "deltaESingle", "deltaLSingle", "deltaJSingle", "deltaMSingle", "deltaNSingle", "deltaEPair", "deltaLPair", "deltaJPair", "deltaMPair", "deltaNPair",
                                   "samebasis", "steps", "precision", "missingCalc", "missingWhittaker", "exponent",
                                   "minEx", "minEy", "minEz", "minBx", "minBy", "minBz", "maxEx", "maxEy", "maxEz", "maxBx", "maxBy", "maxBz", "diamagnetism", "minR", "maxR",
                                   "invE","invO","perE","perO","refE","refO", "conserveM", "sametrafo"]

    # field map of atom 1 and atom 2 (samebasis == True)
    keys_for_cprogram_field12 = ["species1", "n1", "l1", "j1", "m1", "species2", "n2", "l2", "j2", "m2",
                                 "deltaESingle", "deltaLSingle", "deltaJSingle", "deltaMSingle", "deltaNSingle",
                                 "samebasis", "steps", "precision", "missingCalc", "missingWhittaker",
                                 "minEx", "minEy", "minEz", "minBx", "minBy", "minBz", "maxEx", "maxEy", "maxEz", "maxBx", "maxBy", "maxBz", "diamagnetism"]


class MainWindow(QtGui.QMainWindow):

    def __init__(self, parent=None):

        if os.name == 'nt':
            ext = ".exe"
            locale.setlocale(locale.LC_ALL, 'English')
        else:
            ext = ""
            locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

        QtCore.QLocale.setDefault(QtCore.QLocale(locale.getlocale()[0]))

        super().__init__(parent)

        del pg.graphicsItems.GradientEditorItem.Gradients['greyclip']
        del pg.graphicsItems.GradientEditorItem.Gradients['cyclic']
        del pg.graphicsItems.GradientEditorItem.Gradients['spectrum']
        del pg.graphicsItems.GradientEditorItem.Gradients['bipolar']

        # TODO make class, check if module exists
        """from palettable import cubehelix"""

        """color = cubehelix.Cubehelix.make(sat=1.8,n=7,rotation=1.21,start=1.2,reverse=True).colors[::-1]
        color = np.append(color, [[255]]*len(color), axis=1).astype(np.ubyte)"""
        color = [
            [0, 0, 0, 255],
            [12, 67, 0, 255],
            [0, 123, 118, 255],
            [122, 109, 240, 255],
            [255, 121, 197, 255],
            [247, 204, 159, 255],
            [255, 255, 255, 255]
        ]
        pos = np.linspace(0, 1, len(color))
        pg.graphicsItems.GradientEditorItem.Gradients['cubehelix1'] = {
            'mode': 'rgb', 'ticks': [(p, tuple(c)) for c, p in zip(color, pos)]}

        """color = cubehelix.Cubehelix.make(sat=1.5,n=7,rotation=-1.0,start=0.9,reverse=True).colors[::-1]
        color = np.append(color, [[255]]*len(color), axis=1).astype(np.ubyte)"""
        color = [
            [0, 0, 0, 255],
            [75, 19, 77, 255],
            [63, 80, 167, 255],
            [44, 164, 156, 255],
            [117, 206, 113, 255],
            [226, 215, 161, 255],
            [255, 255, 255, 255]
        ]
        pos = np.linspace(0, 1, len(color))
        pg.graphicsItems.GradientEditorItem.Gradients['cubehelix2'] = {
            'mode': 'rgb', 'ticks': [(p, tuple(c)) for c, p in zip(color, pos)]}

        for k, v in pg.graphicsItems.GradientEditorItem.Gradients.items():
            mapticks = [(1 - p, c) for p, c in v['ticks']]
            # mappos = [1-p for p, c in v['ticks']]
            # mapticks[np.argmin(mappos)] = (mapticks[np.argmin(mappos)][0],(245,245,245,255))
            pg.graphicsItems.GradientEditorItem.Gradients[k]['ticks'] = mapticks

        self.ui = Ui_plotwindow()
        self.ui.setupUi(self)

        self.invalidQuantumnumbers = [False, False, False, False]

        self.samebasis = False

        self.systemdict = SystemDict(self.ui)
        self.plotdict = PlotDict(self.ui)

        self.userpath = os.path.expanduser('~')

        self.filepath = self.userpath  # os.getcwd()
        self.systemfile = None
        self.plotfile = None
        self.resultfile = None

        # http://stackoverflow.com/questions/404744/determining-application-path-in-a-python-exe-generated-by-pyinstaller
        if getattr(sys, 'frozen', False):
            self.path_base = os.path.dirname(os.path.realpath(sys.executable))
        elif __file__:
            self.path_base = os.path.dirname(os.path.realpath(__file__))

        if os.path.exists(os.path.join(self.path_base, "conf", "example.sconf")):
            self.path_configurationdir = os.path.join(self.path_base, "conf")
        elif os.path.exists(os.path.join(self.path_base, "../conf", "example.sconf")):
            self.path_configurationdir = os.path.join(self.path_base, "../conf")
        else:
            raise Exception('Directory containing configurations not found.')

        self.path_workingdir = os.path.join(self.path_base, "../../calc")
        self.path_quantumdefects = os.path.join(
            self.path_base, self.path_workingdir, "databases/quantum_defects.db")

        if os.name == 'nt':
            self.path_out = os.path.join(self.userpath, "pairinteraction/")
        else:
            self.path_out = os.path.join(self.userpath, ".pairinteraction/")
        self.path_cache = os.path.join(self.path_out, "cache/")
        self.path_cache_wignerd = os.path.join(self.path_cache, "wignerd/")
        self.path_lastsettings = os.path.join(self.path_out, "lastsettings/")

        self.path_system_last = os.path.join(
            self.path_lastsettings, "lastsettings.sconf")
        self.path_plot_last = os.path.join(
            self.path_lastsettings, "lastsettings.pconf")
        self.path_view_last = os.path.join(
            self.path_lastsettings, "lastview.json")
        self.path_cache_last = os.path.join(
            self.path_lastsettings, "lastcache.json")

        self.path_config = os.path.join(self.path_out, "conf.json")
        self.path_version = os.path.join(self.path_out, "version.json")

        self.proc = None

        self.thread = Worker()

        self.timer = QtCore.QTimer()

        self.momentumcolors = [(55, 126, 184), (77, 175, 74), (228, 26, 28), (152, 78, 163), (
            0, 0, 0), (255 // 5, 255 // 5, 255 // 5)]  # s, p, d, f, other, undetermined

        self.momentummat = [None] * 3
        self.labelmat = [None] * 3
        self.labelstates = [None] * 3
        self.momentumstrings = [None] * 3
        self.stateidx_field = [None] * 3
        self.yMin_field = [None] * 3
        self.yMax_field = [None] * 3

        # fill comboboxes with elements from the quantum defects database
        import sqlite3
        conn = sqlite3.connect(self.path_quantumdefects)
        c = conn.cursor()
        c.execute('SELECT DISTINCT element FROM rydberg_ritz')
        elements = [e[0] for e in c.fetchall()]
        conn.close()

        for combobox in [self.ui.combobox_system_species1, self.ui.combobox_system_species2]:
            combobox.clear()
            combobox.addItems(elements)

        # TODO !!!!!! numBlocks kann auch hoeher als 3 sein!
        self.buffer_basis = [{}, {}, {}]
        self.buffer_energies = [{}, {}, {}]
        self.buffer_positions = [{}, {}, {}]
        self.buffer_boolarr = [{}, {}, {}]
        self.buffer_basis_potential = {}
        self.buffer_energies_potential = {}
        self.buffer_positions_potential = {}

        self.buffer_energiesMap = [{}, {}, {}]
        self.buffer_positionsMap = [{}, {}, {}]
        self.buffer_overlapMap = [{}, {}, {}]
        self.buffer_energiesMap_potential = {}
        self.buffer_positionsMap_potential = {}
        self.buffer_overlapMap_potential = {}

        self.lines_buffer_minIdx = {}
        self.colormap_buffer_minIdx_potential = 0
        self.colormap_buffer_minIdx_field = [0] * 3
        self.lines_buffer_minIdx_field = [0] * 3
        self.iSelected = {}

        """self.ui.colorbutton_plot_nosym.setColor(self.symmetrycolors[0])
        self.ui.colorbutton_plot_sym.setColor(self.symmetrycolors[1])
        self.ui.colorbutton_plot_asym.setColor(self.symmetrycolors[2])"""

        self.ui.colorbutton_plot_invE.setColor((180, 0, 120))
        self.ui.colorbutton_plot_invO.setColor((0, 120, 180))
        self.ui.colorbutton_plot_perE.setColor((180, 60, 60))
        self.ui.colorbutton_plot_perO.setColor((60, 60, 180))
        self.ui.colorbutton_plot_refE.setColor((180, 120, 0))
        self.ui.colorbutton_plot_refO.setColor((120, 0, 180))

        # clrmp = pg.ColorMap(pos,color)
        # self.lut = clrmp.getLookupTable()

        self.ui.gradientwidget_plot_gradient.setOrientation("top")
        self.ui.gradientwidget_plot_gradient.loadPreset('cubehelix1')

        self.tab_field2 = self.ui.tabwidget_plotter.widget(1)

        self.storage_data = [[], [], []]
        self.storage_states = [None, None, None]
        self.storage_configuration = [[None, None], [None, None], [None, None]]

        self.manualRangeX = [False, False, False]
        self.manualRangeY = [False, False, False]

        self.printer = QPrinter(QPrinter.HighResolution)
        self.printer.setPageMargins(20, 15, 20, 20, QPrinter.Millimeter)

        self.momentslabels = ['S', 'P', 'D', 'F']

        self.unperturbedstate = [None, None, None]
        self.overlapstate = [None, None, None]

        self.linesX = [None, None, None]
        self.linesY = [None, None, None]
        self.linesO = [None, None, None]
        self.linesSelected = [0, 0, 0]  # [None,None,None] waere auch ok
        self.linesData = [[], [], []]  # [None,None,None] waere auch ok
        self.linesSender = [None, None, None]

        # TODOs
        self.ui.checkbox_system_override.setEnabled(False)
        self.ui.action_radial_clear.setEnabled(False)
        # self.ui.spinbox_system_exponent.setMaximum(3)
        self.ui.lineedit_system_precision.hide()
        self.ui.label_system_precision.hide()

        # Group up buttons
        self.overlapgroup = QtGui.QButtonGroup()
        self.overlapgroup.addButton(self.ui.radiobutton_plot_overlapDefined)
        self.overlapgroup.addButton(
            self.ui.radiobutton_plot_overlapUnperturbed)

        self.missinggroup = QtGui.QButtonGroup()
        self.missinggroup.addButton(self.ui.radiobutton_system_missingCalc)
        self.missinggroup.addButton(
            self.ui.radiobutton_system_missingWhittaker)
        # self.missinggroup.addButton(self.ui.radiobutton_system_missingError)

        self.matgroup = QtGui.QButtonGroup()
        self.matgroup.addButton(self.ui.radiobutton_system_matCombined)
        self.matgroup.addButton(self.ui.radiobutton_system_matSeparate)

        self.quantizationgroup = QtGui.QButtonGroup()
        self.quantizationgroup.addButton(
            self.ui.radiobutton_system_quantizationZ)
        self.quantizationgroup.addButton(
            self.ui.radiobutton_system_quantizationInteratomic)

        self.symgroup = QtGui.QButtonGroup()
        self.symgroup.addButton(self.ui.radiobutton_symAuto)
        self.symgroup.addButton(self.ui.radiobutton_symManual)

        # Set validators
        validator_double = DoubleValidator()
        validator_doublenone = DoublenoneValidator()
        validator_doublepositive = DoublepositiveValidator()
        validator_doubledelta = DoubledeltaValidator()

        self.ui.lineedit_system_deltaESingle.setValidator(
            validator_doubledelta)
        self.ui.lineedit_system_deltaEPair.setValidator(validator_doubledelta)
        self.ui.lineedit_system_minEx.setValidator(validator_double)
        self.ui.lineedit_system_minEy.setValidator(validator_double)
        self.ui.lineedit_system_minEz.setValidator(validator_double)
        self.ui.lineedit_system_maxEx.setValidator(validator_double)
        self.ui.lineedit_system_maxEy.setValidator(validator_double)
        self.ui.lineedit_system_maxEz.setValidator(validator_double)
        self.ui.lineedit_system_minBx.setValidator(validator_double)
        self.ui.lineedit_system_minBy.setValidator(validator_double)
        self.ui.lineedit_system_minBz.setValidator(validator_double)
        self.ui.lineedit_system_maxBx.setValidator(validator_double)
        self.ui.lineedit_system_maxBy.setValidator(validator_double)
        self.ui.lineedit_system_maxBz.setValidator(validator_double)
        self.ui.lineedit_system_minR.setValidator(validator_doublepositive)
        self.ui.lineedit_system_maxR.setValidator(validator_doublepositive)
        self.ui.lineedit_system_theta.setValidator(validator_double)
        self.ui.lineedit_system_precision.setValidator(
            validator_doublepositive)
        self.ui.lineedit_field1_minE.setValidator(validator_doublenone)
        self.ui.lineedit_field1_maxE.setValidator(validator_doublenone)
        self.ui.lineedit_field2_minE.setValidator(validator_doublenone)
        self.ui.lineedit_field2_maxE.setValidator(validator_doublenone)
        self.ui.lineedit_potential_minE.setValidator(validator_doublenone)
        self.ui.lineedit_potential_maxE.setValidator(validator_doublenone)

        # Connect signals and slots
        self.thread.criticalsignal.connect(self.showCriticalMessage)

        self.ui.spinbox_system_n1.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_system_n2.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_system_l1.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_system_l2.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_system_j1.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_system_j2.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_system_m1.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_system_m2.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_plot_n1.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_plot_n2.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_plot_l1.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_plot_l2.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_plot_j1.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_plot_j2.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_plot_m1.valueChanged.connect(
            self.validateQuantumnumbers)
        self.ui.spinbox_plot_m2.valueChanged.connect(
            self.validateQuantumnumbers)

        self.ui.spinbox_system_j1.editingFinished.connect(
            self.validateHalfinteger)
        self.ui.spinbox_system_j2.editingFinished.connect(
            self.validateHalfinteger)
        self.ui.spinbox_system_m1.editingFinished.connect(
            self.validateHalfinteger)
        self.ui.spinbox_system_m2.editingFinished.connect(
            self.validateHalfinteger)
        self.ui.spinbox_plot_j1.editingFinished.connect(
            self.validateHalfintegerpositiveOrMinusone)
        self.ui.spinbox_plot_j2.editingFinished.connect(
            self.validateHalfintegerpositiveOrMinusone)
        self.ui.spinbox_plot_m1.editingFinished.connect(
            self.validateHalfintegerOrMinusone)
        self.ui.spinbox_plot_m2.editingFinished.connect(
            self.validateHalfintegerOrMinusone)
        self.ui.spinbox_plot_n1.editingFinished.connect(
            self.validateIntegerpositiveOrMinusone)
        self.ui.spinbox_plot_n2.editingFinished.connect(
            self.validateIntegerpositiveOrMinusone)

        self.ui.combobox_system_species1.currentIndexChanged[
            str].connect(self.forbidSamebasis)
        self.ui.combobox_system_species2.currentIndexChanged[
            str].connect(self.forbidSamebasis)

        self.ui.radiobutton_system_pairbasisDefined.toggled.connect(self.togglePairbasis)
        self.ui.radiobutton_plot_overlapDefined.toggled.connect(self.toggleOverlapstate)
        self.ui.radiobutton_plot_log.toggled.connect(self.toggleYScale)
        self.ui.radiobutton_symManual.toggled.connect(self.toggleSymmetrization)

        self.ui.checkbox_plot_antialiasing.toggled.connect(
            self.toggleAntialiasing)
        self.ui.checkbox_system_samebasis.toggled.connect(self.toggleSamebasis)

        self.ui.spinbox_system_deltaNSingle.valueChanged.connect(
            self.adjustPairlimits)
        self.ui.spinbox_system_deltaLSingle.valueChanged.connect(
            self.adjustPairlimits)
        self.ui.spinbox_system_deltaJSingle.valueChanged.connect(
            self.adjustPairlimits)
        self.ui.spinbox_system_deltaMSingle.valueChanged.connect(
            self.adjustPairlimits)

        self.ui.action_system_open.triggered.connect(self.openSystemConf)
        self.ui.action_system_save.triggered.connect(self.saveSystemConf)
        self.ui.action_plot_open.triggered.connect(self.openPlotConf)
        self.ui.action_plot_save.triggered.connect(self.savePlotConf)
        self.ui.action_quit.triggered.connect(self.close)
        self.ui.action_whatsthis.triggered.connect(
            QtGui.QWhatsThis.enterWhatsThisMode)
        self.ui.action_cache_directory.triggered.connect(
            self.changeCacheDirectory)
        self.ui.action_cache_clear.triggered.connect(self.clearCache)
        self.ui.action_print.triggered.connect(self.print)

        self.ui.action_sconf_reset.triggered.connect(self.resetSConf)
        self.ui.action_pconf_reset.triggered.connect(self.resetPConf)

        self.ui.pushbutton_field1_calc.clicked.connect(self.startCalc)
        self.ui.pushbutton_field2_calc.clicked.connect(self.startCalc)
        self.ui.pushbutton_potential_calc.clicked.connect(self.startCalc)

        self.ui.pushbutton_field1_save.clicked.connect(self.saveResult)
        self.ui.pushbutton_field2_save.clicked.connect(self.saveResult)
        self.ui.pushbutton_potential_save.clicked.connect(self.saveResult)

        self.ui.pushbutton_potential_fit.clicked.connect(self.fitC3C6)

        self.timer.timeout.connect(self.checkForData)

        self.ui.graphicsview_field1_plot.sigXRangeChanged.connect(
            self.detectManualRangeX)
        self.ui.graphicsview_field2_plot.sigXRangeChanged.connect(
            self.detectManualRangeX)
        self.ui.graphicsview_potential_plot.sigXRangeChanged.connect(
            self.detectManualRangeX)
        self.ui.graphicsview_field1_plot.sigYRangeChanged.connect(
            self.detectManualRangeY)
        self.ui.graphicsview_field2_plot.sigYRangeChanged.connect(
            self.detectManualRangeY)
        self.ui.graphicsview_potential_plot.sigYRangeChanged.connect(
            self.detectManualRangeY)

        self.ui.spinbox_system_exponent.valueChanged.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_minBx.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_maxBx.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_minBy.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_maxBy.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_minBz.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_maxBz.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_minEx.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_maxEx.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_minEy.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_maxEy.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_minEz.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_maxEz.editingFinished.connect(self.autosetSymmetrization)
        self.ui.lineedit_system_theta.editingFinished.connect(self.autosetSymmetrization)
        self.ui.spinbox_system_m1.valueChanged.connect(self.autosetSymmetrization)
        self.ui.spinbox_system_m2.valueChanged.connect(self.autosetSymmetrization)
        self.ui.checkbox_system_samebasis.stateChanged.connect(self.autosetSymmetrization)


        # Load cache directory
        if os.path.isfile(self.path_cache_last):
            with open(self.path_cache_last, 'r') as f:
                params = json.load(f)
                self.path_cache = params["cachedir"]

        # Check version

        # Load version
        version_settings_saved = None
        version_cache_saved = None
        if os.path.isfile(self.path_version):
            with open(self.path_version, 'r') as f:
                params = json.load(f)
                version_settings_saved = params["version_settings"]
                version_cache_saved = params["version_cache"]

        # Compare version
        '''if os.path.exists(self.path_out) and version_settings_saved != version_settings and version_cache_saved != version_cache: # Poblem: Cachedirectory muss nicht mehr in self.path_out liegen
            msg = QtGui.QMessageBox()
            msg.setText('A new program version has been installed. Due to major changes, cache and settings have to be cleared. This deletes the directory {}.'.format(self.path_out))
            msg.setIcon(QtGui.QMessageBox.Information);
            msg.addButton(QtGui.QMessageBox.Cancel)
            msg.addButton(QtGui.QMessageBox.Ok)
            msg.setDefaultButton(QtGui.QMessageBox.Ok)
            answer = msg.exec()

            # Delete directory
            if answer == QtGui.QMessageBox.Ok:
                shutil.rmtree(self.path_out)
            else:
                sys.exit()'''

        if os.path.exists(self.path_lastsettings) and version_settings_saved != version_settings:
            msg = QtGui.QMessageBox()
            msg.setText('A new program version has been installed. Due to configuration changes, settings have to be cleared. This deletes the directory {}.'.format(
                self.path_lastsettings))
            msg.setIcon(QtGui.QMessageBox.Information)
            msg.addButton(QtGui.QMessageBox.Cancel)
            msg.addButton(QtGui.QMessageBox.Ok)
            msg.setDefaultButton(QtGui.QMessageBox.Ok)
            answer = msg.exec()

            # Delete directory
            if answer == QtGui.QMessageBox.Ok:
                shutil.rmtree(self.path_lastsettings)
            else:
                sys.exit()

        if os.path.exists(self.path_cache) and version_cache_saved != version_cache:
            msg = QtGui.QMessageBox()
            msg.setText('A new program version has been installed. Due to database changes, the cache has to be cleared. This deletes the directory {}.'.format(
                self.path_cache))
            msg.setIcon(QtGui.QMessageBox.Information)
            msg.addButton(QtGui.QMessageBox.Cancel)
            msg.addButton(QtGui.QMessageBox.Ok)
            msg.setDefaultButton(QtGui.QMessageBox.Ok)
            answer = msg.exec()

            # Delete directory
            if answer == QtGui.QMessageBox.Ok:
                shutil.rmtree(self.path_cache)
            else:
                sys.exit()

        # Create directories
        if not os.path.exists(self.path_out):
            os.makedirs(self.path_out)
            if os.name == 'nt':
                import ctypes
                FILE_ATTRIBUTE_HIDDEN = 0x02
                ret = ctypes.windll.kernel32.SetFileAttributesW(
                    self.path_out, FILE_ATTRIBUTE_HIDDEN)
                if not ret:
                    raise ctypes.WinError()

        if not os.path.isfile(self.path_version):
            with open(self.path_version, 'w') as f:
                json.dump({'version_settings': version_settings,
                           'version_cache': version_cache}, f, indent=4, sort_keys=True)

        if not os.path.exists(self.path_lastsettings):
            os.makedirs(self.path_lastsettings)

            with open(self.path_version, 'r') as f:
                version_cache_saved = json.load(f)["version_cache"]

            with open(self.path_version, 'w') as f:
                json.dump({'version_settings': version_settings,
                           'version_cache': version_cache_saved}, f, indent=4, sort_keys=True)

        if not os.path.isfile(self.path_cache_last):
            with open(self.path_cache_last, 'w') as f:
                json.dump({"cachedir": self.path_cache},
                          f, indent=4, sort_keys=True)

        if not os.path.exists(self.path_cache):
            os.makedirs(self.path_cache)

            with open(self.path_version, 'r') as f:
                version_settings_saved = json.load(f)["version_settings"]

            with open(self.path_version, 'w') as f:
                json.dump({'version_settings': version_settings_saved,
                           'version_cache': version_cache}, f, indent=4, sort_keys=True)

        if not os.path.exists(self.path_cache_wignerd):
            os.makedirs(self.path_cache_wignerd)

        # create object to calculate wigner d matrix
        self.wignerd = Wignerd(self.path_cache_wignerd)

        # print(self.wignerd.calc(1/2, 1/2, -1/2, -np.pi/5))

        # Load last settings
        if not os.path.isfile(self.path_system_last):
            shutil.copyfile(os.path.join(
                self.path_configurationdir, "example.sconf"), self.path_system_last)
        self.loadSettingsSystem(self.path_system_last)

        if not os.path.isfile(self.path_plot_last):
            shutil.copyfile(os.path.join(
                self.path_configurationdir, "example.pconf"), self.path_plot_last)
        self.loadSettingsPlotter(self.path_plot_last)

        if os.path.isfile(self.path_view_last):
            self.loadSettingsView(self.path_view_last)

        self.samebasis_state = self.ui.checkbox_system_samebasis.checkState()

        # Emit change-signals in order to let the validation run
        self.ui.spinbox_system_n1.valueChanged.emit(
            self.ui.spinbox_system_n1.value())
        self.ui.spinbox_system_n2.valueChanged.emit(
            self.ui.spinbox_system_n2.value())
        self.ui.spinbox_plot_n1.valueChanged.emit(
            self.ui.spinbox_plot_n1.value())
        self.ui.spinbox_plot_n2.valueChanged.emit(
            self.ui.spinbox_plot_n2.value())

        self.ui.spinbox_system_j1.editingFinished.emit()
        self.ui.spinbox_system_j2.editingFinished.emit()
        self.ui.spinbox_system_m1.editingFinished.emit()
        self.ui.spinbox_system_m2.editingFinished.emit()
        self.ui.spinbox_plot_j1.editingFinished.emit()
        self.ui.spinbox_plot_j2.editingFinished.emit()
        self.ui.spinbox_plot_m1.editingFinished.emit()
        self.ui.spinbox_plot_m2.editingFinished.emit()

        self.ui.combobox_system_species1.currentIndexChanged.emit(
            self.ui.combobox_system_species1.currentIndex())

        self.ui.radiobutton_system_pairbasisDefined.toggled.emit(
            self.ui.radiobutton_system_pairbasisDefined.isChecked())
        self.ui.radiobutton_plot_overlapDefined.toggled.emit(
            self.ui.radiobutton_plot_overlapDefined.isChecked())
        self.ui.radiobutton_plot_log.toggled.emit(
            self.ui.radiobutton_plot_log.isChecked())
        self.ui.radiobutton_symManual.toggled.emit(
            self.ui.radiobutton_symManual.isChecked())

        self.ui.checkbox_plot_antialiasing.toggled.emit(
            self.ui.checkbox_plot_antialiasing.isChecked())
        self.ui.checkbox_system_samebasis.toggled.emit(
            self.ui.checkbox_system_samebasis.isChecked())

        self.ui.spinbox_system_deltaNSingle.valueChanged.emit(
            self.ui.spinbox_system_deltaNSingle.value())
        self.ui.spinbox_system_deltaLSingle.valueChanged.emit(
            self.ui.spinbox_system_deltaLSingle.value())
        self.ui.spinbox_system_deltaJSingle.valueChanged.emit(
            self.ui.spinbox_system_deltaJSingle.value())
        self.ui.spinbox_system_deltaMSingle.valueChanged.emit(
            self.ui.spinbox_system_deltaMSingle.value())

        # Setup plot
        constDistance = self.getConstDistance()
        constEField = self.getConstEField()
        constBField = self.getConstBField()

        self.graphicviews_plot = [self.ui.graphicsview_field1_plot,
                                  self.ui.graphicsview_field2_plot, self.ui.graphicsview_potential_plot]

        for idx in range(3):
            self.graphicviews_plot[idx].setDownsampling(
                ds=True, auto=True, mode='peak')
            self.graphicviews_plot[idx].setClipToView(True)
            self.graphicviews_plot[idx].setLabel(
                'left', 'Energy (' + str(Units.energy) + ')')
            self.graphicviews_plot[idx].scene().contextMenu = None
            self.graphicviews_plot[idx].plotItem.ctrlMenu = None

            self.graphicviews_plot[idx].getAxis("bottom").setZValue(
                1000)  # HACK to bring axis into the foreground
            self.graphicviews_plot[idx].getAxis("left").setZValue(
                1000)  # HACK to bring axis into the foreground

            if (idx in [0, 1] and constEField and not constBField) or (idx == 2 and constDistance and not constBField):
                self.graphicviews_plot[idx].setLabel(
                    'bottom', 'Magnetic field (' + str(Units.bfield) + ')')
            elif (idx in [0, 1]) or (idx == 2 and constDistance and not constEField):
                self.graphicviews_plot[idx].setLabel(
                    'bottom', 'Electric field (' + str(Units.efield) + ')')
            elif (idx == 2):
                self.graphicviews_plot[idx].setLabel(
                    'bottom', 'Interatomic distance (' + str(Units.length) + ')')

    def loadSettingsSystem(self, path):
        with open(path, 'r') as f:
            params = json.load(f)
            for k, v in params.items():
                self.systemdict[k] = Quantity(v[0], v[1])

    def loadSettingsPlotter(self, path):
        with open(path, 'r') as f:
            params = json.load(f)
            self.ui.gradientwidget_plot_gradient.restoreState(
                params["gradientwidget"])
            del params["gradientwidget"]
            for k, v in params.items():
                self.plotdict[k] = Quantity(v[0], v[1])

    def loadSettingsView(self, path):
        with open(path, 'r') as f:
            params = json.load(f)
            self.ui.tabwidget_config.setCurrentIndex(params["config"])
            self.ui.tabwidget_plotter.setCurrentIndex(params["plotter"])
            self.ui.toolbox_system.setCurrentIndex(params["system"])
            if "filepath" in params.keys():
                self.filepath = params["filepath"]

    def saveSettingsSystem(self, path):
        def save(f):
            params = dict()
            for k, v in self.systemdict.items():
                params[k] = [v.magnitude, v.units]
            json.dump(params, f, indent=4, sort_keys=True)

        if isinstance(path, str):
            with open(path, 'w') as f:
                save(f)
        else:
            save(path)

    def saveSettingsPlotter(self, path):
        def save(f):
            params = dict()
            for k, v in self.plotdict.items():
                params[k] = [v.magnitude, v.units]
            params["gradientwidget"] = self.ui.gradientwidget_plot_gradient.saveState()
            json.dump(params, f, indent=4, sort_keys=True)

        if isinstance(path, str):
            with open(path, 'w') as f:
                save(f)
        else:
            save(path)

    def saveSettingsView(self, path):
        def save(f):
            params = dict()
            params["config"] = self.ui.tabwidget_config.currentIndex()
            params["plotter"] = self.ui.tabwidget_plotter.currentIndex()
            params["system"] = self.ui.toolbox_system.currentIndex()
            if self.filepath != self.userpath:
                params["filepath"] = self.filepath
            json.dump(params, f, indent=4, sort_keys=True)

        if isinstance(path, str):
            with open(path, 'w') as f:
                save(f)
        else:
            save(path)

    def get1DPosition(self, vec):
        vec = np.array(vec)
        return np.sign(np.vdot(vec, [1, 1, 1])) * np.linalg.norm(vec)

    def getConstEField(self):
        minVec = np.array([self.systemdict['minEx'].magnitude, self.systemdict[
                          'minEy'].magnitude, self.systemdict['minEz'].magnitude])
        maxVec = np.array([self.systemdict['maxEx'].magnitude, self.systemdict[
                          'maxEy'].magnitude, self.systemdict['maxEz'].magnitude])
        return np.all(minVec == maxVec)

    def getConstBField(self):
        minVec = np.array([self.systemdict['minBx'].magnitude, self.systemdict[
                          'minBy'].magnitude, self.systemdict['minBz'].magnitude])
        maxVec = np.array([self.systemdict['maxBx'].magnitude, self.systemdict[
                          'maxBy'].magnitude, self.systemdict['maxBz'].magnitude])
        return np.all(minVec == maxVec)

    def getConstDistance(self):
        minR = self.systemdict['minR'].magnitude  # TODO
        maxR = self.systemdict['maxR'].magnitude  # TODO
        return minR == maxR

    # def getSameSpecies(self):
    #    return self.systemdict['species1'] == self.systemdict['species2']

    def abortCalculation(self):
        # kill c++ process - this terminates the self.thread, too
        if self.proc is not None:
            self.proc.terminate()

        # wait until self.thread has finished
        self.thread.wait()

        # clear queues
        self.thread.clear()

    def checkForData(self):
        dataamount = 0

        # === print status ===
        elapsedtime = "{}".format(
            timedelta(seconds=int(time() - self.starttime)))
        if self.thread.message != "":
            self.ui.statusbar.showMessage(
                self.thread.message + ", elapsed time " + elapsedtime)
        else:
            self.ui.statusbar.showMessage("Elapsed time " + elapsedtime)

        # === check if memory consumption is to high ===
        if psutil.virtual_memory().percent > 99:  # TODO: is the virtual or swap memory the problem on rqo-donkey?
            self.abortCalculation()
            QtGui.QMessageBox.critical(
                self, "Message", "The program has run out of memory.")

        # === process field and potential maps ===

        for idx in range(3):

            if idx > 0 and not self.thread.dataqueue_field1.empty():
                continue
            if idx > 1 and not self.thread.dataqueue_field2.empty():
                continue

            basisfile = [self.thread.basisfile_field1,
                         self.thread.basisfile_field2, self.thread.basisfile_potential][idx]
            dataqueue = [self.thread.dataqueue_field1,
                         self.thread.dataqueue_field2, self.thread.dataqueue_potential][idx]

            # --- load basis states ---

            if basisfile != "":
                # load basis
                basis = np.loadtxt(basisfile)

                if len(basis.shape) == 1:
                    basis = np.array([basis])

                if basis.size == 0:
                    nState = 0
                else:
                    nState = len(basis)

                if nState == 0:
                    # save basis
                    self.storage_states[idx] = None

                else:
                    # save basis
                    self.storage_states[idx] = basis

                    # determine which state to highlite
                    if self.ui.groupbox_plot_overlap.isChecked():
                        # update status bar
                        message_old = self.ui.statusbar.currentMessage()
                        if idx == 0 and self.thread.samebasis:
                            idxtype = 3
                        else:
                            idxtype = idx
                        status_type = ["Field map of first atom: ", "Field map of second atom: ",
                                       "Pair potential: ", "Field maps: "][idxtype]
                        self.ui.statusbar.showMessage(
                            status_type + "calculate overlap states")
                        QtGui.QApplication.processEvents()

                        # calculate overlap states
                        if self.angle != 0:  # TODO Vereinheitlichen: fuer die verschidenden idx selbe Funktion verwenden, erste Spalte aus basis entfernen
                            if idx == 0:
                                boolarr = self.overlapstate[idx][[0, 1, 2]] != -1
                                stateidx = np.where(np.all(basis[:, [1, 2, 3]][:, boolarr] == self.overlapstate[
                                                    idx][None, [0, 1, 2]][:, boolarr], axis=-1))[0]
                                relevantBasis = basis[stateidx]

                                statecoeff = np.ones_like(stateidx, dtype=np.float)
                                m1 = self.overlapstate[idx][3]
                                if m1 != -1:
                                    for j in np.unique(relevantBasis[:, 3]):
                                        boolarr = np.all(
                                            relevantBasis[:, [3]] == [j], axis=-1)
                                        if np.abs(m1) > j:
                                            statecoeff[boolarr] *= 0
                                        else:
                                            withjBasis = relevantBasis[boolarr]
                                            for m2 in np.unique(withjBasis[:, 4]):
                                                boolarr = np.all(
                                                    relevantBasis[:, [3, 4]] == [j, m2], axis=-1)
                                                statecoeff[
                                                    boolarr] *= self.wignerd.calc(j, m1, m2, self.angle)

                                boolarr = self.overlapstate[
                                    idx][[0, 1, 2, 3]] == -1
                                if sum(boolarr) > 0:
                                    undeterminedQuantumNumbers = relevantBasis[
                                        :, [1, 2, 3, 4]][:, boolarr]
                                    sorter = np.lexsort(
                                        undeterminedQuantumNumbers.T[::-1])
                                    diff = np.append([False], np.diff(
                                        undeterminedQuantumNumbers[sorter], axis=0).any(axis=1))
                                    stateamount = np.cumsum(
                                        diff)[np.argsort(sorter)]
                                else:
                                    stateamount = np.zeros_like(stateidx)

                                if self.thread.samebasis and np.any(self.overlapstate[idx][[0, 1, 2, 3]] != self.overlapstate[idx][[4, 5, 6, 7]]):
                                    boolarr = self.overlapstate[
                                        idx][[4, 5, 6]] != -1
                                    stateidx2 = np.where(np.all(basis[:, [1, 2, 3]][:, boolarr] == self.overlapstate[
                                                         idx][None, [4, 5, 6]][:, boolarr], axis=-1))[0]
                                    relevantBasis = basis[stateidx2]

                                    statecoeff2 = np.ones_like(
                                        stateidx2, dtype=np.float)
                                    m1 = self.overlapstate[idx][7]
                                    if m1 != -1:
                                        for j in np.unique(relevantBasis[:, 3]):
                                            boolarr = np.all(
                                                relevantBasis[:, [3]] == [j], axis=-1)
                                            if np.abs(m1) > j:
                                                statecoeff2[boolarr] *= 0
                                            else:
                                                withjBasis = relevantBasis[boolarr]
                                                for m2 in np.unique(withjBasis[:, 4]):
                                                    boolarr = np.all(
                                                        relevantBasis[:, [3, 4]] == [j, m2], axis=-1)
                                                    statecoeff2[
                                                        boolarr] *= self.wignerd.calc(j, m1, m2, self.angle)

                                    boolarr = self.overlapstate[
                                        idx][[4, 5, 6, 7]] == -1
                                    if sum(boolarr) > 0:
                                        undeterminedQuantumNumbers = relevantBasis[
                                            :, [1, 2, 3, 4]][:, boolarr]
                                        sorter = np.lexsort(
                                            undeterminedQuantumNumbers.T[::-1])
                                        diff = np.append([False], np.diff(
                                            undeterminedQuantumNumbers[sorter], axis=0).any(axis=1))
                                        stateamount2 = np.cumsum(
                                            diff)[np.argsort(sorter)]
                                    else:
                                        stateamount2 = np.zeros_like(stateidx2)

                                    statecoeff = np.append(statecoeff, statecoeff2)
                                    stateidx = np.append(stateidx, stateidx2)
                                    stateamount = np.append(
                                        stateamount, stateamount2)

                                """stateidx = np.where(np.all(basis[:,[1,2,3]] == self.overlapstate[idx][None,[0,1,2]],axis=-1))[0]
                                statecoeff = []
                                j = self.overlapstate[idx][2]
                                for state in basis[stateidx]:
                                    m2 = state[4]
                                    m1 = self.overlapstate[idx][3]
                                    coeff = self.wignerd.calc(j, m2, m1, self.angle)
                                    statecoeff.append(coeff)
                                stateamount = np.zeros_like(stateidx)
                                if self.thread.samebasis and np.any(self.overlapstate[idx][[0,1,2,3]] != self.overlapstate[idx][[4,5,6,7]]):
                                    stateidx_second = np.where(np.all(basis[:,[1,2,3]] == self.overlapstate[idx][None,[4,5,6]],axis=-1))[0]
                                    j = self.overlapstate[idx][6]
                                    for state in basis[stateidx_second]:
                                        m2 = state[4]
                                        m1 = self.overlapstate[idx][7]
                                        coeff = self.wignerd.calc(j, m2, m1, self.angle)
                                        statecoeff.append(coeff)
                                    stateidx = np.append(stateidx, stateidx_second)
                                    stateamount = np.append(stateamount,np.ones_like(stateidx_second))"""
                            elif idx == 1:
                                boolarr = self.overlapstate[idx][[4, 5, 6]] != -1
                                stateidx = np.where(np.all(basis[:, [1, 2, 3]][:, boolarr] == self.overlapstate[
                                                    idx][None, [4, 5, 6]][:, boolarr], axis=-1))[0]
                                relevantBasis = basis[stateidx]

                                statecoeff = np.ones_like(stateidx, dtype=np.float)
                                m1 = self.overlapstate[idx][7]
                                if m1 != -1:
                                    for j in np.unique(relevantBasis[:, 3]):
                                        boolarr = np.all(
                                            relevantBasis[:, [3]] == [j], axis=-1)
                                        if np.abs(m1) > j:
                                            statecoeff[boolarr] *= 0
                                        else:
                                            withjBasis = relevantBasis[boolarr]
                                            for m2 in np.unique(withjBasis[:, 4]):
                                                boolarr = np.all(
                                                    relevantBasis[:, [3, 4]] == [j, m2], axis=-1)
                                                statecoeff[
                                                    boolarr] *= self.wignerd.calc(j, m1, m2, self.angle)

                                boolarr = self.overlapstate[
                                    idx][[4, 5, 6, 7]] == -1
                                if sum(boolarr) > 0:
                                    undeterminedQuantumNumbers = relevantBasis[
                                        :, [1, 2, 3, 4]][:, boolarr]
                                    sorter = np.lexsort(
                                        undeterminedQuantumNumbers.T[::-1])
                                    diff = np.append([False], np.diff(
                                        undeterminedQuantumNumbers[sorter], axis=0).any(axis=1))
                                    stateamount = np.cumsum(
                                        diff)[np.argsort(sorter)]
                                else:
                                    stateamount = np.zeros_like(stateidx)

                                """stateidx = np.where(np.all(basis[:,[1,2,3]] == self.overlapstate[idx][None,[4,5,6]],axis=-1))[0]
                                statecoeff = []
                                j = self.overlapstate[idx][2]
                                for state in basis[stateidx]:
                                    m2 = state[4]
                                    m1 = self.overlapstate[idx][7]
                                    coeff = self.wignerd.calc(j, m2, m1, self.angle)
                                    statecoeff.append(coeff)
                                stateamount = np.zeros_like(stateidx)"""
                            elif idx == 2:
                                boolarr = self.overlapstate[idx][
                                    [0, 1, 2, 4, 5, 6]] != -1
                                stateidx = np.where(np.all(basis[:, [1, 2, 3, 5, 6, 7]][:, boolarr] == self.overlapstate[
                                                    idx][None, [0, 1, 2, 4, 5, 6]][:, boolarr], axis=-1))[0]
                                relevantBasis = basis[stateidx]

                                statecoeff = np.ones_like(stateidx, dtype=np.float)
                                for selector in [0, 4]:
                                    m1 = self.overlapstate[idx][3 + selector]
                                    if m1 != -1:
                                        for j in np.unique(relevantBasis[:, 3 + selector]):
                                            boolarr = np.all(
                                                relevantBasis[:, [3 + selector]] == [j], axis=-1)
                                            if np.abs(m1) > j:
                                                statecoeff[boolarr] *= 0
                                            else:
                                                withjBasis = relevantBasis[boolarr]
                                                for m2 in np.unique(withjBasis[:, 4 + selector]):
                                                    boolarr = np.all(
                                                        relevantBasis[:, [3 + selector, 4 + selector]] == [j, m2], axis=-1)
                                                    statecoeff[
                                                        boolarr] *= self.wignerd.calc(j, m1, m2, self.angle)

                                boolarr = self.overlapstate[idx][
                                    [0, 1, 2, 3, 4, 5, 6, 7]] == -1
                                if sum(boolarr) > 0:
                                    undeterminedQuantumNumbers = relevantBasis[:, [1, 2, 3, 4, 5, 6, 7, 8]][:, boolarr]
                                    sorter = np.lexsort(
                                        undeterminedQuantumNumbers.T[::-1])
                                    diff = np.append([False], np.diff(
                                        undeterminedQuantumNumbers[sorter], axis=0).any(axis=1))
                                    stateamount = np.cumsum(
                                        diff)[np.argsort(sorter)]
                                else:
                                    stateamount = np.zeros_like(stateidx)

                                """stateidx = np.where(np.all(basis[:,[1,2,3,5,6,7]] == self.overlapstate[idx][None,[0,1,2,4,5,6]],axis=-1))[0]
                                statecoeff = []
                                j = self.overlapstate[idx][[2,6]]
                                for state in basis[stateidx]:
                                    m_final = state[[4,8]]
                                    m_initial = self.overlapstate[idx][[3,7]]
                                    coeff = 1
                                    for m2, m1, jj in zip(m_final, m_initial, j):
                                        coeff *= self.wignerd.calc(jj, m2, m1, self.angle)
                                    statecoeff.append(coeff)
                                stateamount = np.zeros_like(stateidx)"""

                            # write calculated wigner d matrix elements into the
                            # cache
                            self.wignerd.save()

                        else:
                            if idx == 0:
                                boolarr = self.overlapstate[idx][[0, 1, 2, 3]] != -1
                                stateidx = np.where(np.all(basis[:, [1, 2, 3, 4]][:, boolarr] == self.overlapstate[idx][None, [0, 1, 2, 3]][:, boolarr], axis=-1))[0]
                                if self.thread.samebasis and np.any(self.overlapstate[idx][[0, 1, 2, 3]] != self.overlapstate[idx][[4, 5, 6, 7]]):
                                    boolarr = self.overlapstate[idx][[4, 5, 6, 7]] != -1
                                    stateidx = np.append(stateidx, np.where(np.all(basis[:, [1, 2, 3, 4]][:, boolarr] == self.overlapstate[idx][None, [4, 5, 6, 7]][:, boolarr], axis=-1))[0])
                            elif idx == 1:
                                boolarr = self.overlapstate[idx][[4, 5, 6, 7]] != -1
                                stateidx = np.where(np.all(basis[:, [1, 2, 3, 4]][:, boolarr] == self.overlapstate[idx][None, [4, 5, 6, 7]][:, boolarr], axis=-1))[0]
                            elif idx == 2:
                                boolarr = self.overlapstate[idx][[0, 1, 2, 3, 4, 5, 6, 7]] != -1
                                stateidx = np.where(np.all(basis[:, [1, 2, 3, 4, 5, 6, 7, 8]][:, boolarr] == self.overlapstate[idx][None, [0, 1, 2, 3, 4, 5, 6, 7]][:, boolarr], axis=-1))[0]

                            statecoeff = np.ones_like(stateidx)
                            stateamount = np.arange(len(stateidx))

                        if len(stateidx) < 1:
                            self.stateidx_field[idx] = None
                        else:
                            self.stateidx_field[idx] = sparse.csc_matrix(
                                (statecoeff, (stateamount, stateidx)), shape=(np.max(stateamount) + 1, nState))

                        # update status bar
                        self.ui.statusbar.showMessage(message_old)
                        QtGui.QApplication.processEvents()

                    else:
                        self.stateidx_field[idx] = None

                    # calculate a matrix that can be used to determine the momenta
                    # inside a basis element
                    if idx == 0 or idx == 1:
                        momentum = basis[:, 2]
                        self.momentummat[idx] = sparse.csc_matrix(
                            (momentum[:, None] == np.arange(np.max(momentum) + 1)[None, :]).astype(int))

                    # extract labels
                    # TODO only if necessary !!!!!

                    if idx == 0 or idx == 1:
                        nlj = basis[:, [1, 2, 3]]
                    elif idx == 2:
                        nlj = basis[:, [1, 2, 3, 5, 6, 7]]

                    # sort pair state names
                    if idx == 2 and self.thread.samebasis:
                        firstsmaller = np.argmax(np.append(nlj[:, 0:3] < nlj[:, 3:6], np.ones(
                            (len(nlj), 1), dtype=bool), axis=-1), axis=-1)  # TODO in Funktion auslagern
                        firstsmaller_reverted = np.argmax(np.append(nlj[:, 3:6] < nlj[:, 0:3], np.ones(
                            (len(nlj), 1), dtype=bool), axis=-1), axis=-1)  # TODO in Funktion auslagern
                        namesToSwap = firstsmaller > firstsmaller_reverted
                        nlj[namesToSwap] = nlj[namesToSwap][:, [3, 4, 5, 0, 1, 2]]

                    sorter = np.lexsort(nlj.T[::-1])
                    nlj = nlj[sorter]
                    diff = np.append([True], np.diff(nlj, axis=0).any(axis=1))
                    cumsum = np.cumsum(diff)[np.argsort(sorter)]

                    # determine labels of states
                    self.labelstates[idx] = nlj[diff]
                    self.labelmat[idx] = sparse.coo_matrix((np.ones_like(cumsum), (np.arange(len(
                        cumsum)), cumsum - 1)), shape=(len(cumsum), len(self.labelstates[idx]))).tocsr()  # nStates, nLabels

                    # determine labels of momenta
                    if idx == 0 or idx == 1:  # TODO !!!
                        self.momentumstrings[idx] = [" {}".format(i) for i in np.arange(
                            np.max(self.labelstates[idx][:, 1]) + 1).astype(np.int)]
                    elif idx == 2:
                        self.momentumstrings[idx] = [" {}".format(i) for i in np.arange(
                            np.max(self.labelstates[idx][:, [1, 4]]) + 1).astype(np.int)]
                    self.momentumstrings[idx][:4] = self.momentslabels

                # remove basis file from hard disk
                os.remove(basisfile)

                # clear variables
                self.lines_buffer_minIdx_field = {}
                self.buffer_basis = {}
                self.buffer_energies = {}
                self.buffer_positions = {}
                self.buffer_boolarr = {}

                self.labelprob = None
                self.labelprob_energy = None

                self.yMin_field[idx] = None
                self.yMax_field[idx] = None

                # indicate that the basis file is already processed
                if idx == 0:
                    self.thread.basisfile_field1 = ""
                elif idx == 1:
                    self.thread.basisfile_field2 = ""
                elif idx == 2:
                    self.thread.basisfile_potential = ""
                    self.thread.basisfile_potential = ""

            # --- check if there is some new data and if yes, plot it ---

            if not dataqueue.empty() and self.storage_states[idx] is None:
                dataqueue.get()  # TODO make this hack unnecessary

            elif not dataqueue.empty():

                graphicsview_plot = [self.ui.graphicsview_field1_plot,
                                     self.ui.graphicsview_field2_plot, self.ui.graphicsview_potential_plot]

                minE = self.minE[idx]
                maxE = self.maxE[idx]

                # --- storage that allows to draw the hole buffer at once, at least if it is no very large ---
                x = np.array([])
                y = np.array([])
                l = np.array([])
                s = []

                while not dataqueue.empty() and dataamount < 5000:  # stop loop if enough data is collected

                    # --- load eigenvalues (energies, y value) and eigenvectors (basis) ---
                    filestep, numBlocks, blocknumber, filename = dataqueue.get()

                    # save data
                    self.storage_data[idx].append(
                        [filestep, blocknumber, filename])

                    eigensystem = Eigensystem(filename)
                    energies = eigensystem.energies
                    basis = eigensystem.basis

                    if idx == 2:
                        symmetrycolor = []

                        if eigensystem.params["inversion"] == "1":
                            symmetrycolor.append(self.ui.colorbutton_plot_invE.color().getRgb()[:-1])
                        elif eigensystem.params["inversion"] == "-1":
                            symmetrycolor.append(self.ui.colorbutton_plot_invO.color().getRgb()[:-1])

                        if eigensystem.params["permutation"] == "1":
                            symmetrycolor.append(self.ui.colorbutton_plot_perE.color().getRgb()[:-1])
                        elif eigensystem.params["permutation"] == "-1":
                            symmetrycolor.append(self.ui.colorbutton_plot_perO.color().getRgb()[:-1])

                        if eigensystem.params["reflection"] == "1":
                            symmetrycolor.append(self.ui.colorbutton_plot_refE.color().getRgb()[:-1])
                        elif eigensystem.params["reflection"] == "-1":
                            symmetrycolor.append(self.ui.colorbutton_plot_refO.color().getRgb()[:-1])

                        if len(symmetrycolor) > 0:
                            symmetrycolor = tuple(np.mean(symmetrycolor, axis=0).astype(int))
                        else:
                            symmetrycolor = (40, 40, 40)

                    # --- determine which basis elements are within the energy range ---
                    boolarr = np.ones(len(energies), dtype=np.bool)
                    boolarr[np.isnan(energies)] = False
                    energies[np.isnan(energies)] = 0

                    if minE is not None:
                        boolarr &= energies >= minE / self.converter_y
                    if maxE is not None:
                        boolarr &= energies <= maxE / self.converter_y

                    # cut the energies
                    energies = energies[boolarr]

                    # convert the energies
                    energies *= self.converter_y

                    # cut the basis
                    idxarr, = np.nonzero(boolarr)
                    transformator = sparse.coo_matrix((np.ones_like(idxarr), (idxarr, np.arange(
                        len(idxarr)))), shape=(basis.shape[1], len(idxarr))).tocsr()
                    basis = basis * transformator

                    # probability to be in a certain state
                    probs = np.abs(basis)  # nState, nBasis
                    probs.data **= 2

                    # --- calculate the momentum that is associated with a basis element ---
                    if idx == 0 or idx == 1:
                        # calculate which momenta appear in a basis element
                        momentum_probabilty = probs.T * \
                            self.momentummat[idx]  # nBasis, nMomentum

                        # keep information about the momentum that appears with
                        # a probability > 0.5 only
                        momentum_probabilty.data[
                            :] *= momentum_probabilty.data > 0.5
                        momentum_probabilty.eliminate_zeros()
                        momentum_probabilty = sparse.coo_matrix(
                            momentum_probabilty)

                        # store the value of this momentum
                        # -1 means no determinable momentum
                        momentum = - \
                            np.ones(momentum_probabilty.shape[0], dtype=np.int)
                        momentum[
                            momentum_probabilty.row] = momentum_probabilty.col

                    # --- calculate the position (x value) ---
                    if self.xAxis[idx] in ['B', 'E']:
                        rotator = np.array([[np.cos(-self.angle), 0, -np.sin(-self.angle)], [
                                           0, 1, 0], [np.sin(-self.angle), 0, np.cos(-self.angle)]])

                    if self.xAxis[idx] == 'B':
                        fields = [[float(eigensystem.params["Bx"])], [float(
                            eigensystem.params["By"])], [float(eigensystem.params["Bz"])]]
                        fields = np.dot(rotator, fields).flatten()
                        position = self.get1DPosition(
                            fields) * self.converter_x[idx]
                    elif self.xAxis[idx] == 'E':
                        fields = [[float(eigensystem.params["Ex"])], [float(
                            eigensystem.params["Ey"])], [float(eigensystem.params["Ez"])]]
                        fields = np.dot(rotator, fields).flatten()
                        position = self.get1DPosition(
                            fields) * self.converter_x[idx]
                    elif self.xAxis[idx] == 'R':
                        position = float(eigensystem.params[
                                         "R"]) * self.converter_x[idx]

                    # --- draw labels at the beginning of the plotting ---
                    if self.ui.groupbox_plot_labels.isChecked() and filestep == 0:

                        # probability to find a label inside a basis element
                        if not hasattr(self, 'labelprob_energy') or self.labelprob_energy is None:
                            # nBasis, nLabels # TODO !!! tocsr
                            self.labelprob = (
                                probs.T * self.labelmat[idx]).tocsr()
                            self.labelprob_energy = [energies]
                        else:
                            csr_vappend(self.labelprob,
                                        (probs.T * self.labelmat[idx]).tocsr())
                            self.labelprob_energy.append(energies)

                        labelprob_num_potential = len(self.labelprob_energy)

                        if labelprob_num_potential == numBlocks:
                            # total probability to find a label
                            cumprob = self.labelprob.sum(axis=0).getA1()
                            boolarr = cumprob > 0.1

                            # normalize in such a way that the total
                            # probability is always one
                            idxarr, = np.nonzero(boolarr)
                            normalizer = sparse.coo_matrix((1 / cumprob[idxarr], (idxarr, np.arange(
                                len(idxarr)))), shape=(self.labelprob.shape[1], len(idxarr))).tocsr()

                            # store the energetic expectation value of the
                            # labels
                            labelenergies = (
                                (self.labelprob * normalizer).T * np.concatenate(self.labelprob_energy))

                            if len(labelenergies) == 0:
                                continue

                            # store the position of the labels
                            labelposition = position

                            # get size and alpha value of labels
                            size = '{}'.format(
                                max(int(round(self.ui.spinbox_plot_szLabel.value() * 11)), 1))
                            alpha = int(
                                round(self.ui.spinbox_plot_transpLabel.value() * 255))

                            # draw the labels
                            for labelstate, labelenergy in zip(self.labelstates[idx][boolarr], labelenergies):

                                if self.leftSmallerRight[idx]:
                                    anchorX = 0
                                else:
                                    anchorX = 1

                                if (idx == 0 and np.all(labelstate == self.unperturbedstate[idx][[0, 1, 2]])) \
                                        or ((idx == 1 or self.thread.samebasis) and np.all(labelstate == self.unperturbedstate[idx][[4, 5, 6]])) \
                                        or (idx == 2 and np.all(labelstate == self.unperturbedstate[idx][[0, 1, 2, 4, 5, 6]])) \
                                        or ((idx == 2 and self.thread.samebasis) and np.all(labelstate == self.unperturbedstate[idx][[4, 5, 6, 0, 1, 2]])):
                                    color_fill = (255, 192, 203, alpha)
                                    color_border = (255, 182, 193, 255)
                                    zvalue = 16
                                else:
                                    color_fill = (250, 235, 215, alpha)
                                    color_border = (255, 228, 181, 255)
                                    zvalue = 15

                                if idx == 0 or idx == 1:
                                    sn, sl, sj = labelstate
                                    text = pg.TextItem(html='<div style="text-align: center; font-size: ' + size + 'pt;">' +
                                                       '<span style="color: rgba(0,0,0,255);">{}{}<sub style="font-size: '.format(
                                                           int(sn), self.momentumstrings[idx][int(sl)]) +
                                                       size + 'pt;">{}/2</sub></span></div>'.format(int(2 * sj)), anchor=(anchorX, 0.5), fill=color_fill, border=color_border)
                                elif idx == 2:
                                    sn1, sl1, sj1, sn2, sl2, sj2 = labelstate
                                    text = pg.TextItem(html='<div style="text-align: center; font-size: ' + size + 'pt;"><span style="color: rgba(0,0,0,255);">' +
                                                       '{}{}<sub style="font-size: '.format(int(sn1), self.momentumstrings[idx][int(sl1)]) +
                                                       size +
                                                       'pt;">{}/2</sub>'.format(
                                                           int(2 * sj1)) +
                                                       ' {}{}<sub style="font-size: '.format(int(sn2), self.momentumstrings[idx][int(sl2)]) +
                                                       size +
                                                       'pt;">{}/2</sub>'.format(
                                                           int(2 * sj2)) +
                                                       '</span></div>', anchor=(anchorX, 0.5), fill=color_fill, border=color_border)

                                text.setPos(labelposition, labelenergy)
                                text.setZValue(zvalue)
                                graphicsview_plot[idx].addItem(text)

                            posx = labelposition * \
                                np.ones_like(labelenergies) + 1e-12
                            posy = labelenergies + 1e-12
                            curve = PointsItem(np.append(
                                posx, posx - 2e-12), np.append(posy, posy - 2e-12), 0, 0, (255, 255, 255))
                            curve.setZValue(5)
                            graphicsview_plot[idx].addItem(curve)

                            # drawing labels take some time
                            dataamount += 3000

                            self.labelprob = None
                            self.labelprob_energy = None

                    # --- draw color map ---
                    if self.ui.groupbox_plot_overlap.isChecked() and self.steps > 1:
                        # --- get settings ---
                        # get size and alpha value
                        size = self.ui.spinbox_plot_szOverlap.value()
                        alpha = min(
                            self.ui.spinbox_plot_transpOverlap.value(), 0.9999)  # HACK

                        # get resolution
                        res = self.ui.spinbox_plot_resolution.value()

                        # calculate size of color map
                        height_pixelunits = res
                        enlargement = int(
                            max(np.round((height_pixelunits / self.steps - 2) / 2), 0))
                        width_pixelunits = 5 + 4 * enlargement

                        # calculate values to smooth the colormap
                        smootherX = (enlargement * 2 + 2) * 1 / 2
                        smootherY = height_pixelunits * 1 / 150 * size

                        # --- build buffers ---
                        # initialize arrays if first run at a new position
                        # ("filestep")
                        if filestep not in self.buffer_positionsMap[idx].keys():
                            self.buffer_positionsMap[idx][filestep] = position
                            self.buffer_energiesMap[idx][filestep] = []
                            self.buffer_overlapMap[idx][filestep] = []

                        # try to get limits
                        if self.yMax_field[idx] is None and maxE is not None:
                            self.yMax_field[idx] = maxE
                        if self.yMin_field[idx] is None and minE is not None:
                            self.yMin_field[idx] = minE

                        # calculate overlap
                        if self.stateidx_field[idx] is not None:
                            overlap = np.abs(self.stateidx_field[idx].conjugate() * basis)
                            overlap.data **= 2
                            overlap = overlap.sum(axis=0).getA1()

                        # check if limits do not exists
                        if self.yMax_field[idx] is None or self.yMin_field[idx] is None:
                            # append the energies to the arrays
                            self.buffer_energiesMap[idx][
                                filestep].append(energies)

                            # append the overlaps to the arrays
                            if self.stateidx_field[idx] is not None:
                                self.buffer_overlapMap[idx][
                                    filestep].append(overlap)
                            else:
                                self.buffer_overlapMap[idx][
                                    filestep].append(np.zeros_like(energies))

                            # check if all data of the zeroth position is
                            # collected
                            # TODO ensure that this also works if numBlocks
                            # changes with the step
                            if 0 in self.buffer_overlapMap[idx].keys() and len(self.buffer_overlapMap[idx][0]) == numBlocks:
                                # make limits
                                if self.yMin_field[idx] is None:
                                    self.yMin_field[idx] = np.nanmin(
                                        np.concatenate(self.buffer_energiesMap[idx][0]))
                                if self.yMax_field[idx] is None:
                                    self.yMax_field[idx] = np.nanmax(
                                        np.concatenate(self.buffer_energiesMap[idx][0]))

                                # calculate energy-indices
                                for f in self.buffer_energiesMap[idx].keys():
                                    for i in range(len(self.buffer_energiesMap[idx][f])):
                                        boolarr = (self.buffer_energiesMap[idx][f][i] >= self.yMin_field[idx]) & (
                                            self.buffer_energiesMap[idx][f][i] <= self.yMax_field[idx])
                                        self.buffer_energiesMap[idx][f][i] = bytescale(self.buffer_energiesMap[idx][f][i][
                                                                                       boolarr], low=0, high=height_pixelunits - 1, cmin=self.yMin_field[idx], cmax=self.yMax_field[idx])
                                        self.buffer_overlapMap[idx][f][
                                            i] = self.buffer_overlapMap[idx][f][i][boolarr]
                        else:
                            # determine whether the energies lie within the
                            # limits
                            boolarr = (energies >= self.yMin_field[idx]) & (
                                energies <= self.yMax_field[idx])

                            # append the energy-indices to the arrays
                            self.buffer_energiesMap[idx][filestep].append(bytescale(
                                energies[boolarr], low=0, high=height_pixelunits - 1, cmin=self.yMin_field[idx], cmax=self.yMax_field[idx]))

                            # append the overlaps to the arrays
                            if self.stateidx_field[idx] is not None:
                                self.buffer_overlapMap[idx][
                                    filestep].append(overlap[boolarr])
                            else:
                                self.buffer_overlapMap[idx][filestep].append(
                                    np.zeros_like(energies[boolarr]))

                        # extract line to fit C3 or C6
                        if self.stateidx_field[idx] is not None and len(overlap) > 0 and (blocknumber not in self.linesO[idx].keys() or np.max(overlap) > 0.5 * np.max(self.linesO[idx][blocknumber])):
                            idxSelected = np.argmax(overlap)
                            xSelected = position
                            eSelected = energies[idxSelected]
                            oSelected = overlap[idxSelected]
                            if blocknumber not in self.linesX[idx].keys():
                                self.linesX[idx][blocknumber] = []
                            if blocknumber not in self.linesY[idx].keys():
                                self.linesY[idx][blocknumber] = []
                            if blocknumber not in self.linesO[idx].keys():
                                self.linesO[idx][blocknumber] = []
                            self.linesX[idx][blocknumber].append(xSelected)
                            self.linesY[idx][blocknumber].append(eSelected)
                            self.linesO[idx][blocknumber].append(oSelected)

                            # the c3 and c6 buttons should be enabled
                            if filestep == 0 and idx == 2:
                                self.ui.pushbutton_potential_fit.setEnabled(
                                    True)
                                self.ui.combobox_potential_fct.setEnabled(True)

                        # --- build color maps starting at the lowest position---
                        # loop over positions ("filestep") as long as three
                        # subsequent positions
                        # ("self.colormap_buffer_minIdx_field[idx]+0,1,2") are
                        # within the buffers
                        while True:
                            # break if limits do not exist
                            if self.yMax_field[idx] is None or self.yMin_field[idx] is None:
                                break

                            # break if buffer index is not in the buffer
                            bufferidx = []

                            # not start
                            if self.colormap_buffer_minIdx_field[idx] != 0:
                                if self.colormap_buffer_minIdx_field[idx] - 1 not in self.buffer_positionsMap[idx].keys():
                                    break
                                bufferidx.append(-1)

                            if self.colormap_buffer_minIdx_field[idx] not in self.buffer_positionsMap[idx].keys():
                                break
                            bufferidx.append(0)

                            # not end
                            if self.colormap_buffer_minIdx_field[idx] != self.steps - 1:
                                if self.colormap_buffer_minIdx_field[idx] + 1 not in self.buffer_positionsMap[idx].keys():
                                    break
                                bufferidx.append(1)

                            # break if the data is not buffered of all blocks,
                            # yet
                            tobreak = False
                            for i in bufferidx:
                                if len(self.buffer_energiesMap[idx][self.colormap_buffer_minIdx_field[idx] + i]) < numBlocks:
                                    tobreak = True
                            if tobreak:
                                break

                            # calculate position-indices
                            positions = [self.buffer_positionsMap[idx][
                                self.colormap_buffer_minIdx_field[idx] + i] for i in bufferidx]

                            # add outer points if at the end or start
                            # end
                            if self.colormap_buffer_minIdx_field[idx] == self.steps - 1:
                                positions = positions + \
                                    [2 * positions[1] - positions[0]]
                            # start
                            elif self.colormap_buffer_minIdx_field[idx] == 0:
                                positions = [2 * positions[0] -
                                             positions[1]] + positions

                            # determine limits of the color map part
                            posLeft = (positions[0] + positions[1]) / 2
                            posRight = (positions[1] + positions[2]) / 2
                            idx_left, idx_right = bytescale(np.array(
                                [posLeft, posRight]), low=0, high=width_pixelunits - 1, cmin=positions[0], cmax=positions[-1])

                            # calculate unit converters
                            self.displayunits2pixelunits_x = width_pixelunits / \
                                (positions[2] - positions[0])
                            self.displayunits2pixelunits_y = height_pixelunits / \
                                (self.yMax_field[idx] - self.yMin_field[idx])

                            # build map
                            # x-y-coordinate system, origin is at the bottom
                            # left corner
                            colormap = np.zeros(
                                (width_pixelunits, height_pixelunits))

                            for i in bufferidx:
                                overlap = np.concatenate(self.buffer_overlapMap[idx][
                                                         self.colormap_buffer_minIdx_field[idx] + i])
                                pos = self.buffer_positionsMap[idx][
                                    self.colormap_buffer_minIdx_field[idx] + i]
                                idx_pos = bytescale(
                                    pos, low=0, high=width_pixelunits - 1, cmin=positions[0], cmax=positions[-1])
                                idx_energies = np.concatenate(self.buffer_energiesMap[idx][
                                                              self.colormap_buffer_minIdx_field[idx] + i])

                                if self.logscale:
                                    overlap[overlap < 1e-2] = 1e-2
                                    overlap = (2 + np.log10(overlap)) / 2

                                colormap[idx_pos, :] = sparse.coo_matrix((overlap, (idx_energies, np.arange(len(
                                    idx_energies)))), shape=(height_pixelunits, len(idx_energies))).sum(axis=1).getA1()

                                dataamount += len(idx_energies) * 10

                            # smoothing
                            colormap = gaussian_filter(
                                colormap, (smootherX, smootherY), mode='constant')

                            # cutting
                            colormap = colormap[idx_left:idx_right]

                            # normalizing
                            normalizer = np.zeros(
                                (width_pixelunits, int(2 * 3 * smootherY + 1)))
                            normalizer[
                                int((width_pixelunits - 1) / 2), int(3 * smootherY)] = 1
                            normalizer = gaussian_filter(
                                normalizer, (smootherX, smootherY), mode='constant')

                            colormap /= np.max(normalizer)

                            # plotting
                            img = pg.ImageItem(
                                image=colormap, opacity=alpha, autoDownsample=True, lut=self.lut, levels=[-0.002, 1])  # HACK
                            img.setRect(QtCore.QRectF(posLeft - 0.5 / self.displayunits2pixelunits_x, self.yMin_field[idx] - 0.5 / self.displayunits2pixelunits_y, posRight - posLeft, self.yMax_field[
                                        idx] - self.yMin_field[idx] + 1 / self.displayunits2pixelunits_y))  # xMin, yMin_field[idx], xSize, ySize # TODO energyMin anpassen wegen Pixelgroesse
                            img.setZValue(3)
                            graphicsview_plot[idx].addItem(img)

                            # remove plotted data from buffer
                            # not start
                            if self.colormap_buffer_minIdx_field[idx] != 0:
                                del self.buffer_energiesMap[idx][
                                    self.colormap_buffer_minIdx_field[idx] - 1]
                                del self.buffer_positionsMap[idx][
                                    self.colormap_buffer_minIdx_field[idx] - 1]
                                del self.buffer_overlapMap[idx][
                                    self.colormap_buffer_minIdx_field[idx] - 1]

                            # increase the buffer index
                            self.colormap_buffer_minIdx_field[idx] += 1

                    # --- draw lines ---
                    if self.ui.groupbox_plot_lines.isChecked():
                        if blocknumber not in self.buffer_basis.keys():
                            self.buffer_basis[blocknumber] = {}
                            self.buffer_energies[blocknumber] = {}
                            self.buffer_positions[blocknumber] = {}
                            self.buffer_boolarr[blocknumber] = {}
                            self.lines_buffer_minIdx_field[blocknumber] = 0
                            """self.iSelected[blocknumber] = None
                            self.linesX[idx][blocknumber] = []
                            self.linesY[idx][blocknumber] = []"""

                        self.buffer_basis[blocknumber][filestep] = basis
                        self.buffer_energies[blocknumber][filestep] = energies
                        self.buffer_positions[blocknumber][filestep] = position
                        self.buffer_boolarr[blocknumber][filestep] = []

                        if idx == 0 or idx == 1:
                            # cut the momenta to reasonable values
                            momentum[momentum > len(
                                self.momentumcolors) - 2] = len(self.momentumcolors) - 2
                            momentum[momentum < 0] = len(
                                self.momentumcolors) - 1

                            # loop over momenta
                            for i in range(len(self.momentumcolors)):
                                # determine which basis elements have the
                                # current momentum
                                boolarr = momentum == i
                                self.buffer_boolarr[blocknumber][
                                    filestep].append(boolarr)
                        elif idx == 2:
                            self.buffer_boolarr[blocknumber][filestep].append(np.ones_like(
                                self.buffer_energies[blocknumber][filestep], dtype=np.bool))

                        # get size and alpha value of points
                        size = self.ui.spinbox_plot_szLine.value()
                        alpha = self.ui.spinbox_plot_transpLine.value() * 255

                        """# legend
                        if idx == 2 and filestep == 0 and symmetry != 0:
                            graphicsview_plot[idx].addLegend()
                            style = pg.PlotDataItem(pen = pg.mkPen(self.symmetrycolors[1]+(alpha,),width=size,cosmetic=True))
                            graphicsview_plot[idx].plotItem.legend.addItem(style, "sym")
                            style = pg.PlotDataItem(pen = pg.mkPen(self.symmetrycolors[2]+(alpha,),width=size,cosmetic=True))
                            graphicsview_plot[idx].plotItem.legend.addItem(style, "asym")"""

                        while self.lines_buffer_minIdx_field[blocknumber] in self.buffer_basis[blocknumber].keys() and self.lines_buffer_minIdx_field[blocknumber] + 1 in self.buffer_basis[blocknumber].keys():
                            # determine the data to plot
                            overlap = np.abs(self.buffer_basis[blocknumber][self.lines_buffer_minIdx_field[blocknumber]].conj(
                            ).T * self.buffer_basis[blocknumber][self.lines_buffer_minIdx_field[blocknumber] + 1])  # nBasis first, nBasis second

                            overlap.data[overlap.data <= np.sqrt(
                                self.ui.spinbox_plot_connectionthreshold.value())] = 0
                            overlap.eliminate_zeros()
                            csr_keepmax(overlap)

                            overlap = overlap.tocoo()

                            iFirst = overlap.row
                            iSecond = overlap.col

                            ydata = np.transpose([self.buffer_energies[blocknumber][self.lines_buffer_minIdx_field[blocknumber]][
                                                 iFirst], self.buffer_energies[blocknumber][self.lines_buffer_minIdx_field[blocknumber] + 1][iSecond]])
                            xdata = np.ones_like(ydata)
                            xdata[:, 0] *= self.buffer_positions[blocknumber][
                                self.lines_buffer_minIdx_field[blocknumber]]
                            xdata[:, 1] *= self.buffer_positions[blocknumber][
                                self.lines_buffer_minIdx_field[blocknumber] + 1]

                            """# track lines with the largest probability to find the overlapstate
                            if self.lines_buffer_minIdx_field[blocknumber] == 0:
                                overlapWithOverlapstate = np.abs(self.stateidx_field[idx].conjugate()*self.buffer_basis[blocknumber][self.lines_buffer_minIdx_field[blocknumber]])
                                overlapWithOverlapstate.data **= 2
                                overlapWithOverlapstate = overlapWithOverlapstate.sum(axis=0).getA1()

                                if len(overlapWithOverlapstate) > 0:
                                    self.iSelected[blocknumber] = np.argmax(overlapWithOverlapstate)
                                    xSelected = self.buffer_positions[blocknumber][self.lines_buffer_minIdx_field[blocknumber]]
                                    eSelected = self.buffer_energies[blocknumber][self.lines_buffer_minIdx_field[blocknumber]][self.iSelected[blocknumber]]
                                    self.linesX[idx][blocknumber].append(xSelected)
                                    self.linesY[idx][blocknumber].append(eSelected)
                                else:
                                    self.iSelected[blocknumber] = -1

                            boolarr = iFirst == self.iSelected[blocknumber]
                            if np.sum(boolarr) == 1:
                                self.iSelected[blocknumber] = iSecond[boolarr]
                                xSelected = self.buffer_positions[blocknumber][self.lines_buffer_minIdx_field[blocknumber]+1]
                                eSelected, = self.buffer_energies[blocknumber][self.lines_buffer_minIdx_field[blocknumber]+1][self.iSelected[blocknumber]]
                                self.linesX[idx][blocknumber].append(xSelected)
                                self.linesY[idx][blocknumber].append(eSelected)
                            else:
                                self.iSelected[blocknumber] = -1"""

                            # loop over momenta
                            numMomenta = len(self.buffer_boolarr[blocknumber][
                                             self.lines_buffer_minIdx_field[blocknumber]])

                            for i in range(numMomenta):
                                boolarr = self.buffer_boolarr[blocknumber][
                                    self.lines_buffer_minIdx_field[blocknumber]][i][iFirst]
                                if np.sum(boolarr) == 0:
                                    continue

                                # determine the associated color
                                if idx == 0 or idx == 1:
                                    color = self.momentumcolors[i]
                                elif idx == 2:
                                    color = symmetrycolor

                                # plot the data
                                # TODO alpha and color der Funktion zusammen
                                # uebergeben
                                curve = MultiLine(xdata[boolarr], ydata[
                                                  boolarr], size, alpha, color)
                                curve.setZValue(7)
                                graphicsview_plot[idx].addItem(curve)

                            del self.buffer_basis[blocknumber][
                                self.lines_buffer_minIdx_field[blocknumber]]
                            del self.buffer_energies[blocknumber][
                                self.lines_buffer_minIdx_field[blocknumber]]
                            del self.buffer_positions[blocknumber][
                                self.lines_buffer_minIdx_field[blocknumber]]
                            del self.buffer_boolarr[blocknumber][
                                self.lines_buffer_minIdx_field[blocknumber]]

                            # increase the buffer index
                            self.lines_buffer_minIdx_field[blocknumber] += 1

                            dataamount += len(iFirst) * 10

                    # --- store data to plot several points at once ---
                    if self.ui.groupbox_plot_points.isChecked():
                        x = np.append(x, position * np.ones_like(energies))
                        y = np.append(y, energies)
                        if idx == 0 or idx == 1:
                            l = np.append(l, momentum)
                        elif idx == 2:
                            s += [symmetrycolor]*len(energies)

                        dataamount += len(x)

                # --- plot the stored data ---
                if self.ui.groupbox_plot_points.isChecked() and len(x) > 0:

                    # get size and alpha value of points
                    size = self.ui.spinbox_plot_szPoint.value()
                    alpha = self.ui.spinbox_plot_transpPoint.value() * 255

                    if idx == 0 or idx == 1:
                        # cut the momenta to reasonable values
                        l[l > len(self.momentumcolors) -
                          2] = len(self.momentumcolors) - 2
                        l[l < 0] = len(self.momentumcolors) - 1

                    # find unique symmetry colors
                    if idx == 0 or idx == 1:
                        looprange = len(self.momentumcolors)
                    elif idx == 2:
                        s = np.array(s)
                        uniquesymmetrycolors = np.unique(s.view(np.dtype((np.void, s.dtype.itemsize*s. shape[1])))).view(s.dtype).reshape(-1, s.shape[1])
                        looprange = len(uniquesymmetrycolors)

                    # loop over momenta
                    for i in range(looprange):
                        if idx == 0 or idx == 1:
                            # determine which basis elements have the current
                            # momentum
                            boolarr = l == i
                            if (np.sum(boolarr) == 0):
                                continue

                            # determine the associated color
                            color = self.momentumcolors[i]
                        elif idx == 2:
                            # determine which basis elements have the current
                            # symmetry
                            boolarr = np.all(s == uniquesymmetrycolors[i], axis=1)
                            if (np.sum(boolarr) == 0):
                                continue

                            # determine the associated color
                            color = tuple(uniquesymmetrycolors[i])

                        # plot the basis elements
                        curve = PointsItem(
                            x[boolarr], y[boolarr], size, alpha, color)
                        curve.setZValue(5)
                        graphicsview_plot[idx].addItem(curve)

                # --- update the graphics view ---
                graphicsview_plot[idx].repaint()

        # check if thread has finished
        if self.thread.isFinished() and self.thread.dataqueue_field1.empty() and self.thread.dataqueue_field2.empty() and self.thread.dataqueue_potential.empty():
            # Delete buffers
            self.buffer_basis = [{}, {}, {}]
            self.buffer_energies = [{}, {}, {}]
            self.buffer_positions = [{}, {}, {}]
            self.buffer_boolarr = [{}, {}, {}]
            self.buffer_basis_potential = {}
            self.buffer_energies_potential = {}
            self.buffer_positions_potential = {}

            self.buffer_energiesMap = [{}, {}, {}]
            self.buffer_positionsMap = [{}, {}, {}]
            self.buffer_overlapMap = [{}, {}, {}]
            self.buffer_energiesMap_potential = {}
            self.buffer_positionsMap_potential = {}
            self.buffer_overlapMap_potential = {}

            self.lines_buffer_minIdx = {}
            self.colormap_buffer_minIdx_potential = 0
            self.colormap_buffer_minIdx_field = [0] * 3
            self.lines_buffer_minIdx_field = [0] * 3

            # Delete c++ process
            if self.proc is not None:
                self.proc.wait()
                self.proc = None

            # Stop this timer
            self.timer.stop()

            # Change buttons
            if self.ui.pushbutton_field1_calc != self.senderbutton:
                self.ui.pushbutton_field1_calc.setEnabled(True)
            if self.ui.pushbutton_field2_calc != self.senderbutton:
                self.ui.pushbutton_field2_calc.setEnabled(True)
            if self.ui.pushbutton_potential_calc != self.senderbutton:
                self.ui.pushbutton_potential_calc.setEnabled(True)
            if self.ui.pushbutton_field1_calc == self.senderbutton:
                self.senderbutton.setText("Calculate field map")
            if self.ui.pushbutton_field2_calc == self.senderbutton:
                self.senderbutton.setText("Calculate field map")
            if self.ui.pushbutton_potential_calc == self.senderbutton:
                self.senderbutton.setText("Calculate potential")

            # Toggle antialiasing # HACK
            if self.ui.checkbox_plot_antialiasing.isChecked():
                for plotarea in [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot, self.ui.graphicsview_potential_plot]:
                    plotarea.setAntialiasing(False)
                    plotarea.setAntialiasing(True)

            # Reset status bar
            self.ui.statusbar.showMessage('')

    @QtCore.pyqtSlot()
    def fitC3C6(self):
        C6notC3 = self.ui.combobox_potential_fct.currentIndex()  # TODO rename variable
        idx = 2
        arrk = list(self.linesX[idx].keys())

        if self.linesSelected[idx] == 0 or self.linesSender[idx] is None or self.linesSender[idx] == C6notC3:
            self.linesSelected[idx] = (
                self.linesSelected[idx] + 1) % (len(arrk) + 1)
        self.linesSender[idx] = C6notC3  # TODO rename variable

        # --- Remove old data from the plot ---
        if len(self.linesData[idx]) > 0:
            for item in self.linesData[idx]:
                self.graphicviews_plot[idx].removeItem(item)
                self.linesData[idx] = []

        # --- Add new data to the plot ---
        if self.linesSelected[idx] > 0:
            k = arrk[self.linesSelected[idx] - 1]

            x = np.array(self.linesX[2][k])
            y = np.array(self.linesY[2][k])

            # HACK: otherwise plot(x,y,...) would not properly work
            sorter = np.argsort(x)
            x = x[sorter]
            y = y[sorter]

            if len(x) >= 2:
                # Get line size
                size = self.ui.spinbox_plot_szLine.value()

                # Plot selected line
                self.linesData[idx].append(self.graphicviews_plot[idx].plot(x, y,
                                                                            pen=pg.mkPen((100, 200, 255, 155), width=size * 2, cosmetic=True)))
                self.linesData[idx][-1].setZValue(13)

                # Plot fitted line
                from scipy import optimize
                y0 = y[np.argmax(x)]
                x0 = np.max(x)
                if C6notC3 == 0:
                    coefftype = ['6']

                    def fitfct(x, c6):
                        return c6 / x**6 - c6 / x0**6 + y0

                elif C6notC3 == 1:
                    coefftype = ['3']

                    def fitfct(x, c3):
                        return c3 / x**3 - c3 / x0**3 + y0

                elif C6notC3 == 2:
                    coefftype = ['3', '6']

                    def fitfct(x, c3, c6):
                        return c3 / x**3 + c6 / x**6 - c3 / x0**3 - c6 / x0**6 + y0

                par, cov = optimize.curve_fit(fitfct, x, y)

                xfit = np.linspace(np.min(x), np.max(x), 300)
                yfit = fitfct(xfit, *par)

                self.linesData[idx].append(self.graphicviews_plot[idx].plot(xfit, yfit,
                                                                            pen=pg.mkPen((0, 200, 200, 255), width=size * 2, style=QtCore.Qt.DotLine, cosmetic=True)))
                self.linesData[idx][-1].setZValue(14)

                # Plot coefficient
                size = '{}'.format(
                    max(int(round(self.ui.spinbox_plot_szLabel.value() * 11)), 1))
                alpha = int(
                    round(self.ui.spinbox_plot_transpLabel.value() * 255))
                color_fill = (200, 255, 255, alpha)
                color_border = (0, 200, 200, 255)

                htmltext = '<div style="text-align: center; font-size: ' + \
                    size + 'pt;"><span style="color: rgba(0,0,0,255);"><b>'
                separator = ''
                for p, c in zip(par, coefftype):
                    htmltext += separator
                    htmltext += '{:.4g} GHz &mu;m<sup style="font-size: '.format(
                        p)
                    htmltext += size + 'pt;">{}</sup>'.format(c)
                    separator = ', '
                htmltext += '</b></span></div>'

                self.linesData[idx].append(pg.TextItem(
                    html=htmltext, fill=color_fill, border=color_border, anchor=(0.5, 0.5)))

                self.linesData[
                    idx][-1].setPos(np.mean(xfit), yfit[np.argmin(np.abs(xfit - np.mean(xfit)))])
                self.linesData[idx][-1].setZValue(20)
                self.graphicviews_plot[idx].addItem(self.linesData[idx][-1])

        # Toggle antialiasing # HACK
        if self.ui.checkbox_plot_antialiasing.isChecked():
            self.graphicviews_plot[idx].setAntialiasing(False)
            self.graphicviews_plot[idx].setAntialiasing(True)

    # @QtCore.pyqtSlot(bool) # TODO !!!!!!!!!!
    def detectManualRangeX(self):
        sender = self.sender()

        idx = -1
        if sender == self.ui.graphicsview_field1_plot:
            idx = 0
        elif sender == self.ui.graphicsview_field2_plot:
            idx = 1
        elif sender == self.ui.graphicsview_potential_plot:
            idx = 2

        if idx > -1:
            self.manualRangeX[idx] = not self.graphicviews_plot[
                idx].getViewBox().getState()["autoRange"][0]

    # @QtCore.pyqtSlot(bool) # TODO !!!!!!!!!!
    def detectManualRangeY(self):
        sender = self.sender()

        idx = -1
        if sender == self.ui.graphicsview_field1_plot.getPlotItem():
            idx = 0
        elif sender == self.ui.graphicsview_field2_plot.getPlotItem():
            idx = 1
        elif sender == self.ui.graphicsview_potential_plot.getPlotItem():
            idx = 2

        if idx > -1:
            self.manualRangeY[idx] = not self.graphicviews_plot[
                idx].getViewBox().getState()["autoRange"][1]

    @QtCore.pyqtSlot(int)
    def adjustPairlimits(self, value):
        sender = self.sender()

        if value == -1:
            maximum = 999
            minimum = -1
        else:
            maximum = value
            minimum = 0

        if sender == self.ui.spinbox_system_deltaNSingle:
            self.ui.spinbox_system_deltaNPair.setMaximum(maximum)
            self.ui.spinbox_system_deltaNPair.setMinimum(minimum)
        elif sender == self.ui.spinbox_system_deltaLSingle:
            self.ui.spinbox_system_deltaLPair.setMaximum(maximum)
            self.ui.spinbox_system_deltaLPair.setMinimum(minimum)
        elif sender == self.ui.spinbox_system_deltaJSingle:
            self.ui.spinbox_system_deltaJPair.setMaximum(maximum)
            self.ui.spinbox_system_deltaJPair.setMinimum(minimum)
        elif sender == self.ui.spinbox_system_deltaMSingle:
            self.ui.spinbox_system_deltaMPair.setMaximum(maximum)
            self.ui.spinbox_system_deltaMPair.setMinimum(minimum)

    @QtCore.pyqtSlot(bool)  # TODO
    def toggleAntialiasing(self):
        checked = self.ui.checkbox_plot_antialiasing.isChecked()
        for plotarea in [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot, self.ui.graphicsview_potential_plot]:
            plotarea.setAntialiasing(checked)
        pg.setConfigOptions(antialias=checked)  # TODO

    @QtCore.pyqtSlot(bool)  # TODO
    def togglePairbasis(self):
        checked = self.ui.radiobutton_system_pairbasisDefined.isChecked()
        self.ui.widget_system_pair.setEnabled(checked)

    @QtCore.pyqtSlot(bool)  # TODO
    def toggleOverlapstate(self):
        checked = self.ui.radiobutton_plot_overlapDefined.isChecked()
        self.ui.widget_plot_qn.setEnabled(checked)

    @QtCore.pyqtSlot(bool)  # TODO
    def toggleSymmetrization(self):
        checked = self.ui.radiobutton_symManual.isChecked()
        self.ui.checkbox_system_invE.setEnabled(checked)
        self.ui.checkbox_system_invO.setEnabled(checked)
        self.ui.checkbox_system_perE.setEnabled(checked)
        self.ui.checkbox_system_perO.setEnabled(checked)
        self.ui.checkbox_system_refE.setEnabled(checked)
        self.ui.checkbox_system_refO.setEnabled(checked)
        self.ui.checkbox_system_conserveM.setEnabled(checked)
        self.autosetSymmetrization()

    def autosetSymmetrization(self):
        if self.ui.radiobutton_symManual.isChecked():
            return

        angle = self.systemdict["theta"].magnitude

        arrlabels = [["minEx", "minEy", "minEz"], ["maxEx", "maxEy", "maxEz"], ["minBx", "minBy", "minBz"], ["maxBx", "maxBy", "maxBz"]]
        rotator = np.array([[np.cos(angle), 0, -np.sin(angle)], [0, 1, 0], [np.sin(angle), 0, np.cos(angle)]])
        fields_unrotated = [np.array([self.systemdict[l].magnitude for l in ls]) for ls in arrlabels]
        fields = [np.dot(rotator, f).flatten() for f in fields_unrotated]

        higherOrder = self.systemdict["exponent"].magnitude > 3
        electricX = fields[0][0] != 0 or fields[1][0] != 0
        electricY = fields[0][1] != 0 or fields[1][1] != 0
        electricZ = fields[0][2] != 0 or fields[1][2] != 0
        magneticX = fields[2][0] != 0 or fields[3][0] != 0
        magneticY = fields[2][1] != 0 or fields[3][1] != 0
        magneticZ = fields[2][2] != 0 or fields[3][2] != 0
        nonzeroM = self.systemdict["m1"].magnitude + self.systemdict["m2"].magnitude != 0
        heteronuclear = not self.systemdict["samebasis"].magnitude

        sym_inversion = (not electricZ) and (not electricX) and (not electricY) and (not heteronuclear)
        sym_permutation = (not higherOrder) and (not heteronuclear)
        sym_reflection = (not magneticZ) and (not magneticX) and (not magneticY) and (not electricX) and (not electricY) and (not nonzeroM)
        sym_rotation = (not magneticX) and (not magneticY) and (not electricX) and (not electricY)

        self.ui.checkbox_system_invE.setChecked(sym_inversion)
        self.ui.checkbox_system_invO.setChecked(sym_inversion)
        self.ui.checkbox_system_perE.setChecked(sym_permutation)
        self.ui.checkbox_system_perO.setChecked(sym_permutation)
        self.ui.checkbox_system_refE.setChecked(sym_reflection)
        self.ui.checkbox_system_refO.setChecked(sym_reflection)
        self.ui.checkbox_system_conserveM.setChecked(sym_rotation)

    @QtCore.pyqtSlot(bool)  # TODO
    def toggleSamebasis(self):
        checked = self.ui.checkbox_system_samebasis.isChecked()
        if checked and self.ui.tabwidget_plotter.count() == 3:
            self.ui.tabwidget_plotter.removeTab(1)
            self.ui.tabwidget_plotter.setTabText(
                0, "Field map of atom 1 and 2")
        elif not checked and self.ui.tabwidget_plotter.count() == 2:
            self.ui.tabwidget_plotter.insertTab(
                1, self.tab_field2, "Field map of atom 2")
            self.ui.tabwidget_plotter.setTabText(0, "Field map of atom 1")

    @QtCore.pyqtSlot(bool)  # TODO
    def toggleYScale(self):
        log = self.ui.radiobutton_plot_log.isChecked()
        if log:
            self.ui.label_plot_1st.setText("< 0.01")  # TODO
            self.ui.label_plot_2nd.setText("0.1")  # TODO
            self.logscale = True
        else:
            self.ui.label_plot_1st.setText("0")
            self.ui.label_plot_2nd.setText("0.5")
            self.logscale = False

    @QtCore.pyqtSlot(str)  # TODO
    def forbidSamebasis(self):
        if self.ui.combobox_system_species1.currentIndex() != self.ui.combobox_system_species2.currentIndex():
            if self.ui.checkbox_system_samebasis.isEnabled():
                self.ui.checkbox_system_samebasis.setEnabled(False)
                self.samebasis_state = self.ui.checkbox_system_samebasis.checkState()
                self.ui.checkbox_system_samebasis.setCheckState(
                    QtCore.Qt.Unchecked)  # TODO !!!!!!!!!!!
        else:
            if not self.ui.checkbox_system_samebasis.isEnabled():
                self.ui.checkbox_system_samebasis.setEnabled(True)
                self.ui.checkbox_system_samebasis.setCheckState(
                    self.samebasis_state)

    @QtCore.pyqtSlot()
    def validateCores(self):
        sender = self.sender()
        if sender.value() != -1 and sender.value() < 2:
            sender.setValue(2)

    @QtCore.pyqtSlot()
    def validateQuantumnumbers(self):
        sender = self.sender()

        quantumnumbers = [
            [self.ui.spinbox_system_n1, self.ui.spinbox_system_l1,
             self.ui.spinbox_system_j1, self.ui.spinbox_system_m1],
            [self.ui.spinbox_system_n2, self.ui.spinbox_system_l2,
             self.ui.spinbox_system_j2, self.ui.spinbox_system_m2],
            [self.ui.spinbox_plot_n1, self.ui.spinbox_plot_l1,
             self.ui.spinbox_plot_j1, self.ui.spinbox_plot_m1],
            [self.ui.spinbox_plot_n2, self.ui.spinbox_plot_l2,
             self.ui.spinbox_plot_j2, self.ui.spinbox_plot_m2]]

        for i, qn in enumerate(quantumnumbers):

            if sender in qn:
                n, l, j, m = qn

                n_err = False
                l_err = False
                j_err = False
                m_err = False

                nn = n.value()
                ll = l.value()
                jj = j.value()
                mm = m.value()

                if (nn != -1 and ll != -1) and (nn - 1 < ll):
                    n_err |= True
                    l_err |= True
                if (ll != -1 and jj != -1) and (abs(ll - jj) != 0.5):
                    l_err |= True
                    j_err |= True
                if (jj != -1 and mm != -1) and (jj < abs(mm)):
                    j_err |= True
                    m_err |= True

                if (nn != -1 and jj != -1) and (nn - 0.5 < jj):
                    n_err |= True
                    j_err |= True
                if (nn != -1 and mm != -1) and (nn - 0.5 < abs(mm)):
                    n_err |= True
                    m_err |= True
                if (ll != -1 and mm != -1) and (ll + 0.5 < abs(mm)):
                    l_err |= True
                    m_err |= True

                if n_err or l_err or j_err or m_err:
                    self.invalidQuantumnumbers[i] = True
                else:
                    self.invalidQuantumnumbers[i] = False

        if np.any(self.invalidQuantumnumbers):
            self.ui.statusbar.showMessage(
                'Invalide quantum numbers specified.')
        else:
            self.ui.statusbar.showMessage('')

    @QtCore.pyqtSlot()
    def validateHalfinteger(self):
        value = self.sender().value()
        self.sender().setValue(np.floor(value) + 0.5)

    @QtCore.pyqtSlot()
    def validateHalfintegerpositiveOrMinusone(self):
        value = self.sender().value()
        if value <= 0:
            self.sender().setValue(-1)
        else:
            self.sender().setValue(np.floor(value) + 0.5)

    @QtCore.pyqtSlot()
    def validateHalfintegerOrMinusone(self):
        value = self.sender().value()
        if value == -1:
            pass
        else:
            self.sender().setValue(np.floor(value) + 0.5)

    @QtCore.pyqtSlot()
    def validateIntegerpositiveOrMinusone(self):
        value = self.sender().value()
        if value <= 0:
            self.sender().setValue(-1)

    @QtCore.pyqtSlot(str)
    def showCriticalMessage(self, msg):
        QtGui.QMessageBox.critical(self, "Message", msg)

    @QtCore.pyqtSlot()
    def startCalc(self):
        if self.proc is None:

            # ensure that validators are called
            focused_widget = QtGui.QApplication.focusWidget()
            if focused_widget is not None:
                focused_widget.clearFocus()

            if np.any(self.invalidQuantumnumbers):
                QtGui.QMessageBox.critical(
                    self, "Message", "Invalide quantum numbers specified.")

            elif self.ui.radiobutton_system_missingWhittaker.isChecked() and max(self.systemdict["n1"].magnitude, self.systemdict["n1"].magnitude) + max(self.systemdict["deltaNSingle"].magnitude, self.systemdict["deltaNPair"].magnitude) > 97:
                QtGui.QMessageBox.critical(
                    self, "Message", "If the principal quantum number exceeds 97, radial matrix elements must be calculated from model potentials.")

            else:
                self.senderbutton = self.sender()

                if self.senderbutton in [self.ui.pushbutton_field1_calc, self.ui.pushbutton_field2_calc] and self.systemdict["theta"].toAU().magnitude != 0:
                    QtGui.QMessageBox.warning(self, "Warning", "For calculating field maps, you might like to set the interaction angle to zero. " +
                                              "A non-zero angle makes the program compute eigenvectors in the rotated basis where the quantization " +
                                              "axis equals the interatomic axis. This slows down calculations.")

                if self.systemdict["theta"].magnitude != 0 and (self.systemdict["deltaMSingle"].magnitude >= 0 or (self.ui.radiobutton_system_pairbasisDefined.isChecked() and self.systemdict["deltaMPair"].magnitude >= 0)):
                    QtGui.QMessageBox.warning(
                        self, "Warning", "For non-zero interaction angles, it is recommended not to restrict the magnetic quantum number.")

                # save last settings
                self.saveSettingsSystem(self.path_system_last)
                self.saveSettingsPlotter(self.path_plot_last)
                self.saveSettingsView(self.path_view_last)

                # change buttons
                if self.senderbutton != self.ui.pushbutton_field1_calc:
                    self.ui.pushbutton_field1_calc.setEnabled(False)
                if self.senderbutton != self.ui.pushbutton_field2_calc:
                    self.ui.pushbutton_field2_calc.setEnabled(False)
                if self.senderbutton != self.ui.pushbutton_potential_calc:
                    self.ui.pushbutton_potential_calc.setEnabled(False)
                if self.senderbutton == self.ui.pushbutton_potential_calc:
                    self.ui.pushbutton_potential_fit.setEnabled(False)
                    self.ui.combobox_potential_fct.setEnabled(False)
                self.senderbutton.setText("Abort calculation")

                # store, whether the same basis should be used for both atoms
                self.samebasis = self.ui.checkbox_system_samebasis.checkState() == QtCore.Qt.Checked

                # store configuration
                # toAU converts the angle from deg to rad
                self.angle = self.systemdict["theta"].toAU().magnitude

                self.minE = [self.plotdict["minE_field1"].magnitude, self.plotdict[
                    "minE_field2"].magnitude, self.plotdict["minE_potential"].magnitude]
                self.maxE = [self.plotdict["maxE_field1"].magnitude, self.plotdict[
                    "maxE_field2"].magnitude, self.plotdict["maxE_potential"].magnitude]

                self.steps = self.systemdict['steps'].magnitude
                # self.calcmissing = self.ui.checkbox_calc_missing.checkState() == QtCore.Qt.Checked

                unperturbedstate = np.array([self.systemdict['n1'].magnitude, self.systemdict['l1'].magnitude, self.systemdict['j1'].magnitude, self.systemdict['m1'].magnitude,
                                             self.systemdict['n2'].magnitude, self.systemdict['l2'].magnitude, self.systemdict['j2'].magnitude, self.systemdict['m2'].magnitude])

                if self.ui.radiobutton_plot_overlapUnperturbed.isChecked():
                    overlapstate = unperturbedstate
                else:
                    overlapstate = np.array([self.plotdict['n1'].magnitude, self.plotdict['l1'].magnitude, self.plotdict['j1'].magnitude, self.plotdict['m1'].magnitude,
                                             self.plotdict['n2'].magnitude, self.plotdict['l2'].magnitude, self.plotdict['j2'].magnitude, self.plotdict['m2'].magnitude])

                self.lut = self.ui.gradientwidget_plot_gradient.getLookupTable(
                    512)

                # clear plots and set them up
                validsenders = [[self.ui.pushbutton_field1_calc, self.ui.pushbutton_potential_calc],
                                [self.ui.pushbutton_field2_calc,
                                    self.ui.pushbutton_potential_calc],
                                [self.ui.pushbutton_potential_calc]]

                constDistance = self.getConstDistance()
                constEField = self.getConstEField()
                constBField = self.getConstBField()

                self.xAxis = [None] * 3
                self.converter_x = [None] * 3
                self.leftSmallerRight = [None] * 3

                for idx in range(3):
                    if (self.senderbutton not in validsenders[idx]) and not (idx == 0 and self.samebasis):
                        continue

                    self.unperturbedstate[idx] = unperturbedstate
                    self.overlapstate[idx] = overlapstate

                    # setup storage variables to save the results
                    filelike_system = StringIO()
                    filelike_plotter = StringIO()
                    self.saveSettingsSystem(filelike_system)
                    self.saveSettingsPlotter(filelike_plotter)

                    self.storage_data[idx] = []
                    self.storage_states[idx] = None
                    self.storage_configuration[idx] = [
                        filelike_system.getvalue(), filelike_plotter.getvalue()]

                    # clear plot
                    autorangestate = self.graphicviews_plot[idx].getViewBox().getState(
                    )["autoRange"]  # HACK to avoid performance issues during clearing
                    if autorangestate[0]:
                        self.graphicviews_plot[idx].disableAutoRange(
                            axis=self.graphicviews_plot[idx].getViewBox().XAxis)
                    if autorangestate[1]:
                        self.graphicviews_plot[idx].disableAutoRange(
                            axis=self.graphicviews_plot[idx].getViewBox().YAxis)
                    self.graphicviews_plot[idx].clear()
                    if autorangestate[0]:
                        self.graphicviews_plot[idx].enableAutoRange(
                            axis=self.graphicviews_plot[idx].getViewBox().XAxis)
                    if autorangestate[1]:
                        self.graphicviews_plot[idx].enableAutoRange(
                            axis=self.graphicviews_plot[idx].getViewBox().YAxis)

                    # set up energy axis
                    self.graphicviews_plot[idx].setLabel(
                        'left', 'Energy (' + str(Units.energy) + ')')
                    self.graphicviews_plot[idx].setLimits(yMin=self.minE[idx])
                    self.graphicviews_plot[idx].setLimits(yMax=self.maxE[idx])

                    # set up step axis
                    if (idx in [0, 1] and constEField and not constBField) or (idx == 2 and constDistance and not constBField):
                        self.xAxis[idx] = 'B'
                        self.graphicviews_plot[idx].setLabel(
                            'bottom', 'Magnetic field (' + str(Units.bfield) + ')')
                        self.converter_x[idx] = Quantity(
                            1, Units.au_bfield).toUU().magnitude
                        posMin = self.get1DPosition([self.systemdict['minBx'].magnitude, self.systemdict[
                                                    'minBy'].magnitude, self.systemdict['minBz'].magnitude])
                        posMax = self.get1DPosition([self.systemdict['maxBx'].magnitude, self.systemdict[
                                                    'maxBy'].magnitude, self.systemdict['maxBz'].magnitude])
                    elif (idx in [0, 1]) or (idx == 2 and constDistance and not constEField):
                        self.xAxis[idx] = 'E'
                        self.graphicviews_plot[idx].setLabel(
                            'bottom', 'Electric field (' + str(Units.efield) + ')')
                        self.converter_x[idx] = Quantity(
                            1, Units.au_efield).toUU().magnitude
                        posMin = self.get1DPosition([self.systemdict['minEx'].magnitude, self.systemdict[
                                                    'minEy'].magnitude, self.systemdict['minEz'].magnitude])
                        posMax = self.get1DPosition([self.systemdict['maxEx'].magnitude, self.systemdict[
                                                    'maxEy'].magnitude, self.systemdict['maxEz'].magnitude])
                    elif (idx == 2):
                        self.xAxis[idx] = 'R'
                        self.graphicviews_plot[idx].setLabel(
                            'bottom', 'Interatomic distance (' + str(Units.length) + ')')
                        self.converter_x[idx] = Quantity(
                            1, Units.au_length).toUU().magnitude
                        posMin = self.systemdict['minR'].magnitude
                        posMax = self.systemdict['maxR'].magnitude

                    self.leftSmallerRight[idx] = posMin < posMax

                    # enable / disable auto range
                    if self.ui.checkbox_plot_autorange.isChecked():
                        self.graphicviews_plot[idx].enableAutoRange()
                    else:
                        if not self.manualRangeX[idx]:
                            self.graphicviews_plot[
                                idx].setXRange(posMin, posMax)
                        if not self.manualRangeY[idx] and self.minE[idx] is not None and self.maxE[idx] is not None:
                            self.graphicviews_plot[idx].setYRange(
                                self.minE[idx], self.maxE[idx])

                    # clear variables
                    self.linesSelected[idx] = 0
                    self.linesData[idx] = []
                    self.linesX[idx] = {}
                    self.linesY[idx] = {}
                    self.linesO[idx] = {}
                    self.linesSender[idx] = None

                self.converter_y = Quantity(
                    1, Units.au_energy).toUU().magnitude

                # save configuration to json file
                with open(self.path_config, 'w') as f:
                    if self.senderbutton == self.ui.pushbutton_potential_calc:
                        keys = self.systemdict.keys_for_cprogram_potential
                    elif self.samebasis:
                        keys = self.systemdict.keys_for_cprogram_field12
                    elif self.senderbutton == self.ui.pushbutton_field1_calc:
                        keys = self.systemdict.keys_for_cprogram_field1
                    elif self.senderbutton == self.ui.pushbutton_field2_calc:
                        keys = self.systemdict.keys_for_cprogram_field2

                    params = {k: self.systemdict[
                        k].toAU().magnitude for k in keys}

                    if self.senderbutton == self.ui.pushbutton_potential_calc:
                        params["zerotheta"] = self.angle == 0

                    if params["deltaNSingle"] < 0:
                        params["deltaNSingle"] = -1
                    if params["deltaLSingle"] < 0:
                        params["deltaLSingle"] = -1
                    if params["deltaJSingle"] < 0:
                        params["deltaJSingle"] = -1
                    if params["deltaMSingle"] < 0:
                        params["deltaMSingle"] = -1

                    if self.senderbutton == self.ui.pushbutton_potential_calc and self.ui.radiobutton_system_pairbasisSame.isChecked():
                        params["deltaNPair"] = -1
                        params["deltaLPair"] = -1
                        params["deltaJPair"] = -1
                        params["deltaMPair"] = -1

                    if self.angle != 0:
                        arrlabels = [["minEx", "minEy", "minEz"], ["maxEx", "maxEy", "maxEz"], [
                            "minBx", "minBy", "minBz"], ["maxBx", "maxBy", "maxBz"]]
                        rotator = np.array([[np.cos(self.angle), 0, -np.sin(self.angle)], [
                                           0, 1, 0], [np.sin(self.angle), 0, np.cos(self.angle)]])

                        for labels in arrlabels:
                            fields = np.array([[params[label]]
                                               for label in labels])
                            fields = np.dot(rotator, fields).flatten()
                            for field, label in zip(fields, labels):
                                params[label] = field

                    """if self.angle != 0 or params["minEx"] != 0 or params["minEy"] != 0 or params["maxEx"] != 0 or params["maxEy"] != 0 or params["minBx"] != 0 or params["minBy"] != 0 or params["maxBx"] != 0 or params["maxBy"] != 0:
                        params["conserveM"] = False
                    else:
                        params["conserveM"] = True

                    if (self.senderbutton == self.ui.pushbutton_potential_calc and params["exponent"] > 3) or params["minEx"] != 0 or params["minEy"] != 0 or params["minEz"] != 0 or params["maxEx"] != 0 or params["maxEy"] != 0 or params["maxEz"] != 0:
                        params["conserveParityL"] = False
                    else:
                        params["conserveParityL"] = True"""

                    # TODO make quantities of None type accessible without
                    # .magnitude

                    # TODO remove this hack
                    if self.senderbutton == self.ui.pushbutton_potential_calc and params["exponent"] == 3:
                        params["dd"] = True
                        params["dq"] = False
                        params["qq"] = False
                    # TODO remove this hack
                    elif self.senderbutton == self.ui.pushbutton_potential_calc and params["exponent"] == 2:
                        params["dd"] = False
                        params["dq"] = False
                        params["qq"] = False
                    else:  # TODO remove this hack
                        params["dd"] = True
                        params["dq"] = True
                        params["qq"] = True

                    json.dump(params, f, indent=4, sort_keys=True)

                # start c++ process
                if params["minEy"] != 0 or params["maxEy"] != 0 or \
                        params["minBy"] != 0 or params["maxBy"] != 0:
                    pi = pic
                else:
                    pi = pir

                self.numprocessors = self.systemdict["cores"].magnitude
                # OMP_NUM_THREADS – Specifies the number of threads to
                # use in parallel regions.  The value of this variable
                # shall be a comma-separated list of positive
                # integers; the value specified the number of threads
                # to use for the corresponding nested level.  If
                # undefined one thread per CPU is used.
                pi.thread_ctrl(self.numprocessors)

                self.starttime = time()


                class Communicator:
                    def __init__(self):
                        self.context = zmq.Context()
                        self.socket = self.context.socket(zmq.SUB)
                        self.socket.bind("tcp://*:5556")
                        self.socket.setsockopt_string(zmq.SUBSCRIBE, u"")

                    def __iter__(self):
                        return self

                    def __next__(self):
                        string = self.socket.recv_string()

                        if ">>END" in string:
                            self.socket.close()
                            self.context.destroy()
                            raise StopIteration
                        else:
                            return string

                # start thread that collects the output
                t = multiprocessing.Process(target=pi.compute,
                                            args=("/home/user/.pairinteraction/conf.json",
                                                  "/home/user/.pairinteraction/cache/"))
                t.start()
                self.thread.execute(Communicator())

                # start timer used for processing the results
                self.timer.start(0)

        else:
            self.abortCalculation()

    @QtCore.pyqtSlot()
    def saveResult(self):
        senderbutton = self.sender()
        if senderbutton == self.ui.pushbutton_field1_save:
            idx = 0
        elif senderbutton == self.ui.pushbutton_field2_save:
            idx = 1
        elif senderbutton == self.ui.pushbutton_potential_save:
            idx = 2

        # get filename
        path = self.plotfile if self.plotfile is not None else self.filepath

        if idx == 0 and self.ui.checkbox_system_samebasis.isChecked():
            description = "field map of atom 1 and 2"
        else:
            description = ["field map of atom 1",
                           "field map of atom 2", "pair potential"][idx]

        filename, _ = QtGui.QFileDialog.getSaveFileName(self,
                                                        "Save {}".format(description), path, "zip (*.zip)")

        if not filename:
            return

        self.resultfile = filename
        self.filepath = os.path.dirname(filename)

        self.saveToZipfile(filename, idx)

    def saveToZipfile(self, filename, idx):
        # open zip file
        ziparchive = zipfile.ZipFile(
            filename, 'w', compression=zipfile.ZIP_STORED)  # zipfile.ZIP_DEFLATED

        # TODO sicherstellen, dass die Funktion auch ohne Plot / mit leerem
        # Plotwindow nicht abst?rzt

        try:
            # save plot
            plotitem = [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot,
                        self.ui.graphicsview_potential_plot][idx].getPlotItem()
            exporter = exporters.ImageExporter(plotitem)
            exporter.parameters()['width'] = 2000
            exporter.parameters()['height'] = 2000
            exporter.parameters()['antialias'] = True
            image = exporter.export(toBytes=True)

            buffer = QtCore.QBuffer()
            buffer.open(QtCore.QIODevice.WriteOnly)
            image = image.scaled(
                1500, 1500, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
            image.save(buffer, "PNG")
            ziparchive.writestr('plot.png', buffer.data())

            # save configuration
            ziparchive.writestr(
                'settings.sconf', self.storage_configuration[idx][0])
            ziparchive.writestr(
                'settings.pconf', self.storage_configuration[idx][1])

            # create data dictionary
            data = {}

            data['numStates'] = 0
            data['numSteps'] = 0
            data['numEigenvectors'] = 0
            data['numOverlapvectors'] = 0
            data['states'] = []
            data['eigenvectors'] = []
            data['eigenvalues'] = []
            data['overlapvectors'] = []
            data['overlaps'] = []
            data['bfields'] = []
            data['efields'] = []

            if idx == 2:
                data['distances'] = []

            if idx == 0 or idx == 1:
                data[
                    'states_description'] = 'state(idxState, {idxState, n, l, j, m})'
            elif idx == 2:
                data[
                    'states_description'] = 'state(idxState, {idxState, n1, l1, j1, m1, n2, l2, j2, m2})'
            data[
                'eigenvectors_description'] = 'eigenvectorcoordinate(0, idxStep)(idxState, idxEigenvector)'
            data[
                'eigenvalues_description'] = 'eigenvalue(idxStep, idxEigenvector)'
            data[
                'overlapvectors_description'] = 'overlapvectorcoordinate(idxOverlapvector, idxState)'
            data['overlaps_description'] = 'overlap(idxStep, idxEigenvector)'
            data['bfields_description'] = 'bfield(idxStep, {Bx, By, Bz})'
            data['efields_description'] = 'efield(idxStep, {Ex, Ey, Ez})'
            if idx == 2:
                data['distances_description'] = 'distance(0, idxStep)'

            # save states
            if self.storage_states[idx] is not None:
                # nState, i-n1-l1-j1-m1-n2-l2-j2-m2 # nState, i-n-l-j-m
                data['states'] = self.storage_states[idx]
                data['numStates'] = len(data['states'])

            # save overlaps
            if self.stateidx_field[idx] is not None:
                data['numOverlapvectors'] = self.stateidx_field[idx].shape[0]
                data['overlapvectors'] = self.stateidx_field[idx]

            # save data
            # TODO Variablen an anderer Stelle anlegen
            self.converter_bfield = \
                Quantity(1, Units.au_bfield).toUU().magnitude
            self.converter_efield = \
                Quantity(1, Units.au_efield).toUU().magnitude
            self.converter_length = \
                Quantity(1, Units.au_length).toUU().magnitude

            rotator = np.array([[np.cos(-self.angle), 0, -np.sin(-self.angle)],
                                [0, 1, 0],
                                [np.sin(-self.angle), 0, np.cos(-self.angle)]])  # TODO !!!!!!!!! self.angle[idx] verwenden

            filestep_last = None

            for filestep, blocknumber, filename in sorted(self.storage_data[idx], key=itemgetter(0, 1)):
                eigensystem = Eigensystem(filename)
                energies = eigensystem.energies * self.converter_y  # nBasis
                # nState, nBasis (stored in Compressed Sparse Column format,
                # CSC)
                basis = eigensystem.basis

                if self.stateidx_field[idx] is not None:
                    overlaps = np.abs(self.stateidx_field[
                                      idx].conjugate() * basis)
                    overlaps.data **= 2
                    overlaps = overlaps.sum(axis=0).getA1()

                if filestep != filestep_last:  # new step
                    field = np.array([float(eigensystem.params["Bx"]), float(
                        eigensystem.params["By"]), float(eigensystem.params["Bz"])])
                    data['bfields'].append(
                        np.dot(rotator, field).flatten() * self.converter_bfield)
                    field = np.array([float(eigensystem.params["Ex"]), float(
                        eigensystem.params["Ey"]), float(eigensystem.params["Ez"])])
                    data['efields'].append(
                        np.dot(rotator, field).flatten() * self.converter_efield)
                    if idx == 2:
                        data['distances'].append(
                            float(eigensystem.params["R"]) * self.converter_length)
                    data['eigenvectors'].append(basis)
                    data['eigenvalues'].append(energies)
                    if self.stateidx_field[idx] is not None:
                        data['overlaps'].append(overlaps)
                else:  # new block
                    csc_happend(data['eigenvectors'][-1], basis)
                    data[
                        'eigenvalues'][-1] = np.append(data['eigenvalues'][-1], energies)
                    if self.stateidx_field[idx] is not None:
                        data[
                            'overlaps'][-1] = np.append(data['overlaps'][-1], overlaps)

                filestep_last = filestep

            if len(data['eigenvalues']) > 0:
                data['numSteps'] = len(data['eigenvalues'])
                data['numEigenvectors'] = len(data['eigenvalues'][0])

            if self.ui.radiobutton_system_quantizationZ.isChecked() and self.angle != 0:

                # find relevant states
                statesum = np.zeros(data['numStates'], dtype=np.float)
                for s in range(data['numSteps']):
                    statesum += np.abs(data['eigenvectors']
                                       [s]).sum(axis=1).getA1()
                statesum += np.abs(data['overlapvectors']).sum(axis=0).getA1()

                boolarr = statesum > 0
                idxconverter = np.arange(data['numStates'])[boolarr]
                relevantstates = data['states'][boolarr]

                # build transformator
                if len(data['states'][0]) == 5:
                    shifts = [0]
                else:
                    shifts = [0, 4]

                for selector in shifts:
                    idx1 = []
                    idx2 = []
                    val = []

                    for j in np.unique(relevantstates[:, 3 + selector]):
                        arrM = np.unique(relevantstates[:, 4 + selector])

                        for m1 in arrM:  # m_row
                            if np.abs(m1) > j:
                                continue
                            boolarr1 = np.all(
                                relevantstates[:, [3 + selector, 4 + selector]] == [j, m1], axis=-1)
                            idxconverter1 = idxconverter[boolarr1]

                            for m2 in arrM:  # m_col
                                if np.abs(m2) > j:
                                    continue
                                boolarr2 = np.all(
                                    relevantstates[:, [3 + selector, 4 + selector]] == [j, m2], axis=-1)
                                idxconverter2 = idxconverter[boolarr2]

                                nl1 = relevantstates[boolarr1][
                                    :, [1 + selector, 2 + selector]]
                                nl2 = relevantstates[boolarr2][
                                    :, [1 + selector, 2 + selector]]
                                matches = sparse.coo_matrix(
                                    np.all(nl1[:, None, :] == nl2[None, :, :], axis=-1))

                                wignerd = self.wignerd.calc(
                                    j, m2, m1, -self.angle)

                                idx1 += idxconverter1[matches.row].tolist()
                                idx2 += idxconverter2[matches.col].tolist()
                                val += [wignerd] * matches.nnz

                    if selector == 0:
                        transformator = sparse.csc_matrix(
                            (val, (idx1, idx2)), shape=(data['numStates'], data['numStates']))
                    else:
                        transformator = transformator.multiply(sparse.csc_matrix(
                            (val, (idx1, idx2)), shape=(data['numStates'], data['numStates'])))

                # write calculated wigner d matrix elements into the cache
                self.wignerd.save()

                # apply transformator
                for s in range(data['numSteps']):
                    data['eigenvectors'][s] = transformator * \
                        data['eigenvectors'][s]

                # print(np.round(np.real(data['eigenvectors'][0][boolarr][:,:5].todense()),2))

            if self.systemdict["matCombined"].magnitude is True:
                filelike = BytesIO()
                io.savemat(filelike, data, do_compression=False,
                           format='5', oned_as='row')
                ziparchive.writestr('data.mat', filelike.getvalue())
            else:
                for idxStep in range(data['numSteps']):
                    data_stepwise = {}
                    data_stepwise['numStates'] = data['numStates']
                    data_stepwise['numSteps'] = 1
                    data_stepwise['numEigenvectors'] = data['numEigenvectors']
                    data_stepwise['numOverlapvectors'] = data[
                        'numOverlapvectors']
                    data_stepwise['states'] = data['states']
                    data_stepwise['eigenvectors'] = [
                        data['eigenvectors'][idxStep]]
                    data_stepwise['eigenvalues'] = [
                        data['eigenvalues'][idxStep]]
                    data_stepwise['overlapvectors'] = data['overlapvectors']
                    data_stepwise['overlaps'] = [data['overlaps'][idxStep]]
                    data_stepwise['bfields'] = [data['bfields'][idxStep]]
                    data_stepwise['efields'] = [data['efields'][idxStep]]

                    if idx == 2:
                        data_stepwise['distances'] = [
                            data['distances'][idxStep]]

                    if idx == 0 or idx == 1:
                        data_stepwise['states_description'] = \
                            'state(idxState, {idxState, n, l, j, m})'

                    elif idx == 2:
                        data_stepwise['states_description'] = \
                            'state(idxState, {idxState, n1, l1, j1, m1, n2, l2, j2, m2})'
                    data_stepwise['eigenvectors_description'] = \
                        'eigenvectorcoordinate(0, 0)(idxState, idxEigenvector)'
                    data_stepwise['eigenvalues_description'] = \
                        'eigenvalue(0, idxEigenvector)'
                    data_stepwise['overlapvectors_description'] = \
                        'overlapvectorcoordinate(idxOverlapvector, idxState)'
                    data_stepwise['overlaps_description'] = \
                        'overlap(0, idxEigenvector)'
                    data_stepwise['bfields_description'] = \
                        'bfield(0, {Bx, By, Bz})'
                    data_stepwise['efields_description'] = \
                        'efield(0, {Ex, Ey, Ez})'

                    if idx == 2:
                        data_stepwise['distances_description'] = \
                            'distance(0, 0)'

                    filelike = BytesIO()
                    io.savemat(filelike, data_stepwise,
                               do_compression=False, format='5', oned_as='row')
                    ziparchive.writestr('data_{:04d}.mat'.format(
                        idxStep), filelike.getvalue())

        finally:
            # close zip file
            ziparchive.close()

        """ symmetry = eigensystem.params["symmetry"] # TODO ausgelesene Symmetrie auch beim Plotten verwenden !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                # save configuration
                if firstround:
                    firstround = False

                    del eigensystem.params["Bx"]
                    del eigensystem.params["By"]
                    del eigensystem.params["Bz"]
                    del eigensystem.params["Ex"]
                    del eigensystem.params["Ey"]
                    del eigensystem.params["Ez"]
                    if idx == 2:
                        del eigensystem.params["R"]
                        del eigensystem.params["symmetry"]

                    # TODO abgespeicherte Konfiguration ladbar machen

                    print(eigensystem.params)"""

    @QtCore.pyqtSlot()
    def saveSystemConf(self):
        path = self.systemfile if self.systemfile is not None else self.filepath
        filename, _ = QtGui.QFileDialog.getSaveFileName(
            self, "Save system configuration", path, "sconf (*.sconf)")

        if filename:
            self.saveSettingsSystem(filename)
            self.systemfile = filename
            self.filepath = os.path.dirname(filename)

    @QtCore.pyqtSlot()
    def savePlotConf(self):
        path = self.plotfile if self.plotfile is not None else self.filepath
        filename, _ = QtGui.QFileDialog.getSaveFileName(
            self, "Save plot configuration", path, "pconf (*.pconf)")

        if filename:
            self.saveSettingsPlotter(filename)
            self.plotfile = filename
            self.filepath = os.path.dirname(filename)

    @QtCore.pyqtSlot()
    def changeCacheDirectory(self):
        text, ok = QtGui.QInputDialog.getText(self, 'Input Dialog',
                                              'Enter new cache directory (the original directory is not deleted):', QtGui.QLineEdit.Normal, self.path_cache)

        if ok:
            self.path_cache = text
            if not os.path.exists(self.path_cache):
                os.makedirs(self.path_cache)

            # Save cache directory
            with open(self.path_cache_last, 'w') as f:
                json.dump({"cachedir": self.path_cache},
                          f, indent=4, sort_keys=True)

    @QtCore.pyqtSlot()
    def clearCache(self):
        files = ['cache_elements.db', 'cache_matrix_complex.db', 'cache_matrix_real.db', 'cache_matrix_complex',
                 'cache_matrix_real']  # TODO: sicherstellen, dass gleiche Namen wie im C++ Programm
        for file in files:
            path = os.path.join(self.path_cache, file)
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)

        if os.path.isdir(self.path_cache_wignerd):
            shutil.rmtree(self.path_cache_wignerd)
        os.makedirs(self.path_cache_wignerd)

    @QtCore.pyqtSlot()
    def openSystemConf(self):
        filename, _ = QtGui.QFileDialog.getOpenFileName(self,
                                                        "Open system configuration", self.filepath, "sconf (*.sconf)")

        if filename:
            self.loadSettingsSystem(filename)
            self.systemfile = filename
            self.filepath = os.path.dirname(filename)

    @QtCore.pyqtSlot()
    def openPlotConf(self):
        filename, _ = QtGui.QFileDialog.getOpenFileName(self,
                                                        "Open plot configuration", self.filepath, "pconf (*.pconf)")

        if filename:
            self.loadSettingsPlotter(filename)
            self.plotfile = filename
            self.filepath = os.path.dirname(filename)

    @QtCore.pyqtSlot()
    def resetSConf(self):
        conf_used = self.path_system_last
        conf_original = os.path.join(self.path_configurationdir, "example.sconf")

        if os.path.isfile(conf_used):
            os.remove(conf_used)
        shutil.copyfile(conf_original, conf_used)
        self.loadSettingsSystem(conf_used)

    @QtCore.pyqtSlot()
    def resetPConf(self):
        conf_used = self.path_plot_last
        conf_original = os.path.join(self.path_configurationdir, "example.pconf")

        if os.path.isfile(conf_used):
            os.remove(conf_used)
        shutil.copyfile(conf_original, conf_used)
        self.loadSettingsPlotter(conf_used)

    @QtCore.pyqtSlot()
    def print(self):
        idx = self.ui.tabwidget_plotter.currentIndex()
        if idx == 1 and self.ui.tabwidget_plotter.count() == 2:
            idx = 2

        # TODO initialize these variables already in the constructor
        if self.storage_configuration[idx][0] is None:
            filelike = StringIO()
            self.saveSettingsSystem(filelike)
            self.storage_configuration[idx][0] = filelike.getvalue()

        if self.storage_configuration[idx][1] is None:
            filelike = StringIO()
            self.saveSettingsPlotter(filelike)
            self.storage_configuration[idx][1] = filelike.getvalue()

        if self.unperturbedstate[idx] is None:
            self.unperturbedstate[idx] = np.array([self.systemdict['n1'].magnitude, self.systemdict['l1'].magnitude, self.systemdict['j1'].magnitude, self.systemdict['m1'].magnitude,
                                                   self.systemdict['n2'].magnitude, self.systemdict['l2'].magnitude, self.systemdict['j2'].magnitude, self.systemdict['m2'].magnitude])

        if self.overlapstate[idx] is None:
            if self.ui.radiobutton_plot_overlapUnperturbed.isChecked():
                self.overlapstate[idx] = self.unperturbedstate[idx]
            else:
                self.overlapstate[idx] = np.array([self.plotdict['n1'].magnitude, self.plotdict['l1'].magnitude, self.plotdict['j1'].magnitude, self.plotdict['m1'].magnitude,
                                                   self.plotdict['n2'].magnitude, self.plotdict['l2'].magnitude, self.plotdict['j2'].magnitude, self.plotdict['m2'].magnitude])

        dialog = QPrintDialog(self.printer, self)
        if dialog.exec_() == QPrintDialog.Accepted:

            tabname = self.ui.tabwidget_plotter.tabText(idx)

            # TODO check if results exist !!!!!!!!!!!!!!!!!!
            # TODO reihenfolge plot settings umdrehen  !!!!!!!!!!!!!!!!!!

            '''if idx == 2: pairpotential = 1
            else: pairpotential = 0'''

            painter = QtGui.QPainter(self.printer)

            # printer settings
            spacer = QtCore.QRectF(painter.viewport()).height() / 40
            font = QtGui.QFont("Helvetica", 10)
            margin = 30
            penwidth = 30

            # load configuration
            sconf = json.loads(self.storage_configuration[idx][0])
            pconf = json.loads(self.storage_configuration[idx][1])

            su = [""] * 8
            st = self.unperturbedstate[idx]
            for k in [0, 4]:
                su[k] = "{}".format(int(st[k]))
            for k in [1, 5]:
                if int(st[k]) >= 0 and int(st[k]) < len(self.momentslabels):
                    su[k] = self.momentslabels[int(st[k])]
                else:
                    su[k] = "{}".format(int(st[k]))
            for k in [2, 3, 6, 7]:
                if float(st[k]).is_integer():
                    ssuo[k] = "{}".format(int(st[k]))
                else:
                    su[k] = "{}/2".format(int(2 * float(st[k])))

            so = [""] * 8
            st = self.overlapstate[idx]
            for k in [0, 4]:
                so[k] = "{}".format(int(st[k]))
            for k in [1, 5]:
                if int(st[k]) >= 0 and int(st[k]) < len(self.momentslabels):
                    so[k] = self.momentslabels[int(st[k])]
                else:
                    so[k] = "{}".format(int(st[k]))
            for k in [2, 3, 6, 7]:
                if float(st[k]).is_integer():
                    so[k] = "{}".format(int(st[k]))
                else:
                    so[k] = "{}/2".format(int(2 * float(st[k])))

            if sconf["pairbasisSame"][0]:
                sconf["deltaNPair"][0] = -1
                sconf["deltaLPair"][0] = -1
                sconf["deltaJPair"][0] = -1
                sconf["deltaMPair"][0] = -1

            interaction = "multipole expansion up to order {}".format(sconf[
                                                                      "exponent"][0])

            # text formating
            '''if idx in [0,2] or sconf["pairbasisSame"][0]:
                u1l = "<b>"
                u1r = "</b>"
            else:
                u1l = ""
                u1r = ""

            if idx in [1,2] or sconf["pairbasisSame"][0]:
                u2l = "<b>"
                u2r = "</b>"
            else:
                u2l = ""
                u2r = ""'''

            u1l = ""
            u1r = ""
            u2l = ""
            u2r = ""

            # image
            plotitem = [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot,
                        self.ui.graphicsview_potential_plot][idx].getPlotItem()
            exporter = exporters.ImageExporter(plotitem)
            exporter.parameters()['antialias'] = True

            # colormap
            arr = self.ui.gradientwidget_plot_gradient.getLookupTable(
                501, alpha=False)
            arr[0] = 0
            arr[250] = 0
            arr[-1] = 0
            bgra = np.empty((1, 501, 4), np.uint8, 'C')
            bgra[..., 0] = arr[..., 2]
            bgra[..., 1] = arr[..., 1]
            bgra[..., 2] = arr[..., 0]
            image_colormap = QtGui.QImage(
                bgra.data, 501, 1, QtGui.QImage.Format_RGB32)
            image_colormap.ndarray = bgra

            # --- generate description ---
            rect = QtCore.QRectF(painter.viewport())
            doc = QtGui.QTextDocument()
            doc.documentLayout().setPaintDevice(self.printer)
            doc.setDefaultFont(font)
            doc.setPageSize(rect.size())

            text = ""

            # text += "<h3>System configuration</h3>"

            # state
            text += "<h4>Unperturbed state</h4>"
            text += "<table cellspacing=5><tr>"
            text += "<td>| " + u1l + "{} {} {}<sub>{}</sub>, m={}".format(sconf["species1"][0], su[0], su[1], su[2], su[3]) + u1r + "; " +\
                u2l + "{} {} {}<sub>{}</sub>, m={}".format(sconf["species2"][0], su[4], su[
                                                           5], su[6], su[7]) + u2r + " &gt;</td>"
            text += "</tr></table>"

            '''text += "unperturbed state "
            text += "| "+u1l+"{} {} {}<sub>{}/2</sub>, m={}/2".format(sconf["species1"][0],sconf["n2"][0],sconf["l1"][0],int(sconf["j1"][0]*2),int(sconf["m1"][0]*2))+u1r+"; "+\
                u2l+"{} {} {}<sub>{}/2</sub>, m={}/2".format(sconf["species2"][0],sconf["n2"][0],sconf["l2"][0],int(sconf["j2"][0]*2),int(sconf["m2"][0]*2))+u2r+" >"'''

            # basis
            text += "<h4>Basis (use same basis for both atoms: {})</h4>".format(
                {True: "yes", False: "no"}[sconf["samebasis"][0]])
            text += "<table cellspacing=5 width=100%><tr>"
            text += "<td>single atom state:</td>"
            text += "<td>&Delta;E = {} {}</td>".format(
                sconf["deltaESingle"][0], sconf["deltaESingle"][1])
            text += "<td>&Delta;n = {}</td>".format(sconf["deltaNSingle"][0])
            text += "<td>&Delta;l = {}</td>".format(sconf["deltaLSingle"][0])
            text += "<td>&Delta;j = {}</td>".format(sconf["deltaJSingle"][0])
            text += "<td>&Delta;m = {}</td>".format(sconf["deltaMSingle"][0])
            if idx == 2:
                text += "</tr><tr>"
            if idx == 2:
                text += "<td>pair state:</td>"
            if idx == 2:
                text += "<td>&Delta;E = {} {}</td>".format(
                    sconf["deltaEPair"][0], sconf["deltaEPair"][1])
            if idx == 2:
                text += "<td>&Delta;n = {}</td>".format(sconf["deltaNPair"][0])
            if idx == 2:
                text += "<td>&Delta;l = {}</td>".format(sconf["deltaLPair"][0])
            if idx == 2:
                text += "<td>&Delta;j = {}</td>".format(sconf["deltaJPair"][0])
            if idx == 2:
                text += "<td>&Delta;m = {}</td>".format(sconf["deltaMPair"][0])
            text += "</tr></table>"

            # interaction
            text += "<h4>Interaction (resolution: {} steps)</h4>".format(
                sconf["steps"][0])
            text += "<table cellspacing=5 width=100%><tr>"
            text += "<td>at start:</td>"
            text += "<td>E = ({}, {}, {}) {}</td>".format(sconf["minEx"][0], sconf[
                "minEy"][0], sconf["minEz"][0], sconf["minEz"][1])
            text += "<td>B = ({}, {}, {}) {}</td>".format(sconf["minBx"][0], sconf[
                "minBy"][0], sconf["minBz"][0], sconf["minBz"][1])
            if idx == 2:
                text += "<td>R = {} {}</td>".format(
                    sconf["minR"][0], sconf["minR"][1])
            text += "</tr><tr>"
            text += "<td>at end:</td>"
            text += "<td>E = ({}, {}, {}) {}</td>".format(sconf["maxEx"][0], sconf[
                "maxEy"][0], sconf["maxEz"][0], sconf["maxEz"][1])
            text += "<td>B = ({}, {}, {}) {}</td>".format(sconf["maxBx"][0], sconf[
                "maxBy"][0], sconf["maxBz"][0], sconf["maxBz"][1])
            if idx == 2:
                text += "<td>R = {} {}</td>".format(
                    sconf["maxR"][0], sconf["maxR"][1])
            text += "</tr></table>"

            if idx == 2:
                text += "<table cellspacing=5><tr>"
            if idx == 2:
                text += "<td>{}, interaction angle {} {}</td>".format(
                    interaction, sconf["theta"][0], sconf["theta"][1])
            if idx == 2:
                text += "</tr></table>"

            # text += "<h3>{}</h3>".format(["Field map","Pair potential"][1])

            doc.setHtml(text)
            height_doc = doc.documentLayout().documentSize().height()

            # --- paint header ---
            painter.save()
            rect = QtCore.QRectF(painter.viewport())
            doc_header = QtGui.QTextDocument()
            doc_header.documentLayout().setPaintDevice(self.printer)
            doc_header.setDefaultFont(font)
            doc_header.setPageSize(rect.size())

            text = ""

            # header
            # text += "<table cellspacing=5 width='100%' bgcolor='#f0f0f0'><tr>"
            # text += "<td><h1>{}</h1></td>".format(["Field map of atom 1","Field map of atom 2","Pair potential","Field map of atom 1 and 2"][idx])
            # text += "</tr></table>"

            text += "<table width=100%><tr>"
            text += "<td align=left>Rydberg interaction calculator</td>"
            text += "<td align=center></td>"
            text += "<td align=right>{}</td>".format(
                datetime.today().strftime("%x"))
            text += "</tr></table>"
            text += "<h2>{}</h2>".format(tabname)

            doc_header.setHtml(text)
            height_header = doc_header.documentLayout().documentSize().height()
            doc_header.drawContents(painter)
            painter.restore()

            # --- paint image ---

            # pen = QtGui.QPen(QtGui.QColor(220, 220, 220), penwidth)
            # pen.setJoinStyle(QtCore.Qt.RoundJoin);
            # painter.setPen(pen)
            # painter.drawRect(rect)

            painter.save()
            rect = painter.viewport()
            rect = QtCore.QRect(rect.x(), rect.y() + height_header + 300, rect.width(
            ), rect.height() - height_doc - height_header - 300 - 400 - spacer * 0.7 - 520)
            exporter.parameters()['width'] = int(rect.width() / 4)
            exporter.parameters()['height'] = int(rect.height() / 4)
            image = exporter.export(toBytes=True)
            size = image.size()
            size.scale(rect.size(), QtCore.Qt.KeepAspectRatio)

            border = QtCore.QRect(rect.x() - 15, rect.y() - 15,
                                  size.width() + 30, size.height() + 30 + spacer * 0.7 + 520)
            pen = QtGui.QPen(QtGui.QColor(220, 220, 220), 70)
            painter.setPen(pen)
            painter.drawRect(border)
            pen = QtGui.QPen(QtGui.QColor(0, 0, 0), 10)
            painter.setPen(pen)
            painter.drawRect(border)
            pen = QtGui.QPen(QtGui.QColor(0, 0, 0), 0)
            painter.setPen(pen)
            brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
            painter.setBrush(brush)
            painter.drawRect(border)

            painter.setViewport(rect.x(), rect.y(),
                                size.width(), size.height())
            painter.setWindow(image.rect())
            painter.drawImage(0, 0, image)
            height_image = size.height()
            width_image = size.width()
            height_box = size.height() + 30 + spacer * 0.7 + 520
            painter.restore()

            # --- paint colormap ---

            painter.save()
            rect = painter.viewport()
            rect = QtCore.QRect(rect.x() + 70, rect.y() + height_image +
                                height_header + 300 + spacer * 0.7 + 220, width_image - 140, 100)
            size = image_colormap.size()
            size.scale(rect.size(), QtCore.Qt.IgnoreAspectRatio)

            painter.save()
            painter.setViewport(rect.x(), rect.y(),
                                size.width(), size.height())
            painter.setWindow(image_colormap.rect())
            painter.drawImage(0, 0, image_colormap)
            painter.restore()

            border = QtCore.QRect(rect.x(), rect.y(),
                                  size.width(), size.height())
            pen = QtGui.QPen(QtGui.QColor(100, 100, 100), 5)
            painter.setPen(pen)
            painter.drawRect(border)
            painter.restore()

            painter.save()
            doc_label = QtGui.QTextDocument()
            doc_label.documentLayout().setPaintDevice(self.printer)
            doc_label.setDefaultFont(font)
            rect = QtCore.QRectF(rect.x(), rect.y() - 220, rect.width(), 999)
            doc_label.setPageSize(rect.size())

            text = ""
            text += "Overlap with | {} {} {}<sub>{}</sub>, m={}; {} {} {}<sub>{}</sub>, m={} &gt;".\
                    format(sconf["species1"][0], so[0], so[1], so[2], so[3],
                           sconf["species2"][0], so[4], so[5], so[6], so[7])
            text += "<br>"
            text += "<table cellspacing=0 width=100%><tr>"
            if pconf["log"][0]:
                text += "<td align=left>&lt; 0.01</td><td align=center>0.1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td><td align=right>1</td>"
            else:
                text += "<td align=left>0</td><td align=center>0.5</td><td align=right>1</td>"
            text += "</tr></table>"

            # text += "<h3>{}</h3>".format(["Field map","Pair potential"][1])

            doc_label.setHtml(text)
            painter.translate(rect.x(), rect.y())
            doc_label.drawContents(painter)
            painter.restore()

            # --- paint description ---
            painter.save()
            painter.translate(0, height_box + 300 + height_header + 400)
            doc.drawContents(painter)
            painter.restore()

            """
            painter.save()
            rect = painter.viewport()
            rect = QtCore.QRect(rect.x(),rect.y()+height_doc+spacer*0.8,rect.width(),height_image+spacer*0.4)
            pen = QtGui.QPen(QtGui.QColor(0, 0, 0), 10)
            painter.setPen(pen)
            painter.drawLine(rect.topLeft(),rect.topRight())
            painter.drawLine(rect.bottomLeft(),rect.bottomRight())
            painter.restore()"""

            """rect = painter.viewport().size()
            size = image.size()

            size.scale(QtCore.QSize(rect.width()-2*margin,rect.height()-2*margin), QtCore.Qt.KeepAspectRatio)
            painter.setViewport(margin, height_doc+spacer+margin, size.width(), size.height())
            painter.setWindow(image.rect())
            painter.drawImage(0, 0, image)

            painter.restore()


            painter.save()


            pen = QtGui.QPen(QtGui.QColor(220, 220, 220), penwidth)
            #pen = QtGui.QPen(QtGui.QColor(0, 0, 0), penwidth)
            pen.setJoinStyle(QtCore.Qt.RoundJoin);
            painter.setPen(pen)


            painter.drawRect(QtCore.QRectF(penwidth/2+5, height_doc+spacer+penwidth/2+5, size.width()+2*margin-penwidth/2-5, size.height()+2*margin-penwidth/2-5))

            #bound = QtCore.QRectF(5, height_doc+spacer+5, size.width()+2*margin-5, size.height()+2*margin-5-penwidth/2)
            #painter.drawLine(bound.bottomLeft(),bound.bottomRight())
            #painter.drawLine(bound.topLeft(),bound.topRight())

            #print(bound.bottomLeft(),bound.bottomRight())
            #print(bound.bottomLeft(),bound.bottomRight())

            plotheight = size.height()+2*margin"""

            """# --- paint settings ---

            painter.save()

            doc = QtGui.QTextDocument()

            rect_doc = QtCore.QRectF(painter.viewport())
            rect_doc = QtCore.QRectF(0,plotheight+height_doc+2*spacer,rect_doc.width(),rect_doc.height()-plotheight-height_doc-2*spacer)
            #painter.drawRect(rect_doc)



            doc.documentLayout().setPaintDevice(self.printer)
            doc.setDefaultFont(font)
            doc.setPageSize(rect_doc.size())


            #ziparchive.writestr('settings.sconf', self.storage_configuration[idx][0])
            #ziparchive.writestr('settings.pconf', self.storage_configuration[idx][1])

            sconf = json.loads(self.storage_configuration[idx][0])
            pconf = json.loads(self.storage_configuration[idx][1])

            s1 = sconf["species1"][0]
            n1 = sconf["n1"][0]
            l1 = sconf["l1"][0]
            j1 = sconf["j1"][0]
            m1 = sconf["m1"][0]
            s2 = sconf["species2"][0]
            n2 = sconf["n2"][0]
            l2 = sconf["l2"][0]
            j2 = sconf["j2"][0]
            m2 = sconf["m2"][0]
            minEx = sconf["minEx"]
            minEy = sconf["minEy"]
            minEz = sconf["minEz"]
            maxEx = sconf["maxEx"]
            maxEy = sconf["maxEy"]
            maxEz = sconf["maxEz"]
            minBx = sconf["minBx"]
            minBy = sconf["minBy"]
            minBz = sconf["minBz"]
            maxBx = sconf["maxBx"]
            maxBy = sconf["maxBy"]
            maxBz = sconf["maxBz"]
            minR = sconf["minR"]
            maxR = sconf["maxR"]
            theta = sconf["theta"]
            steps = sconf["steps"][0]

            if l1 < len(self.momentslabels): l1 = self.momentslabels[l1]
            if l2 < len(self.momentslabels): l2 = self.momentslabels[l2]

            text = ""

            text += "<h4>Unperturbed state</h4>"
            text += "<table cellspacing=5 style='width:100%'><tr>"
            text += "<td>| {} {} {}<sub>{}/2</sub>, m={}/2; {} {} {}<sub>{}/2</sub>, m={}/2 ></td>".format(s1,n1,l1,int(j1*2),int(m1*2),s2,n2,l2,int(j2*2),int(m2*2))
            text += "</tr></table>"

            text += "<h4>Basis</h4>"
            text += "<table cellspacing=5 style='width:100%'><tr>"
            text += "<td>&Delta;E<sub>single</sub> = {} {}</td>".format(minEz[0], minEz[1])
            text += "<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>"
            text += "<td>&Delta;n<sub>single</sub> = {}</td>".format(0)
            text += "<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>"
            text += "<td>&Delta;l<sub>single</sub> = {}</td>".format(0)
            text += "<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>"
            text += "<td>&Delta;j<sub>single</sub> = {}</td>".format(0)
            text += "<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>"
            text += "<td>&Delta;m<sub>single</sub> = {}</td>".format(0)
            text += "</tr><tr>"
            text += "<td>&Delta;E<sub>pair</sub> = {} {}</td>".format(maxEz[0], maxEz[1])
            text += "<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>"
            text += "<td>&Delta;n<sub>pair</sub> = {}</td>".format(0)
            text += "<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>"
            text += "<td>&Delta;l<sub>pair</sub> = {}</td>".format(0)
            text += "<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>"
            text += "<td>&Delta;j<sub>pair</sub> = {}</td>".format(0)
            text += "<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>"
            text += "<td>&Delta;m<sub>pair</sub> = {}</td>".format(0)
            text += "</tr></table>"

            text += "<h4>Fields and interaction</h4>"
            text += "<table cellspacing=5 style='width:100%'><tr>"
            text += "<td>E<sub>start</sub> = ({}, {}, {}) {}</td>".format(minEx[0], minEy[0], minEz[0], minEz[1])
            text += "<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>"
            text += "<td>E<sub>end</sub> = ({}, {}, {}) {}</td>".format(maxEx[0], maxEy[0], maxEz[0], maxEz[1])
            text += "</tr><tr>"
            text += "<td>B<sub>start</sub> = ({}, {}, {}) {}</td>".format(minBx[0], minBy[0], minBz[0], minBz[1])
            text += "<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>"
            text += "<td>B<sub>end</sub> = ({}, {}, {}) {}</td>".format(maxBx[0], maxBy[0], maxBz[0], maxBz[1])
            text += "</tr><tr>"
            text += "<td>R<sub>start</sub> = {} {}</td>".format(minR[0], minR[1])
            text += "<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>"
            text += "<td>R<sub>end</sub> = {} {}</td>".format(maxR[0], maxR[1])
            text += "</tr></table>"

            text += "<table cellspacing=5 style='width:100%'><tr>"
            text += "<td>interaction angle: {} {}</td>".format(theta[0], theta[1])
            text += "</tr><tr>"
            text += "<td>interaction order: dipole-dipole, dipole-quadrupole, quadrupole-quadrupole</td>".format(theta[0], theta[1])
            text += "</tr><tr>"
            text += "<td>steps: {}</td>".format(steps)
            text += "</tr></table>"




            doc.setHtml(text)
            painter.translate(rect_doc.left(), rect_doc.top());
            doc.drawContents(painter)

            painter.restore()"""

            """rect = painter.viewport()
            size = self.image.size()
            size.scale(rect.size(), Qt.KeepAspectRatio)
            painter.setViewport(rect.x(), rect.y(), size.width(), size.height())
            painter.setWindow(self.image.rect())
            painter.drawImage(0, 0, self.image)"""
            """rect = painter.viewport()
            size = self.imageLabel.pixmap().size()
            size.scale(rect.size(), Qt.KeepAspectRatio)
            painter.setViewport(rect.x(), rect.y(), size.width(), size.height())
            painter.setWindow(self.imageLabel.pixmap().rect())
            painter.drawPixmap(0, 0, self.imageLabel.pixmap())
            painter.end()

            # --- paint image ---

            painter.save()

            plotitem = [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot, self.ui.graphicsview_potential_plot][idx].getPlotItem()
            exporter = pg.exporters.ImageExporter(plotitem)
            exporter.parameters()['width'] = 2000
            exporter.parameters()['height'] = 2000
            exporter.parameters()['antialias'] = True
            image = exporter.export(toBytes=True)

            margin = 30
            rect = painter.viewport().size()
            size = image.size()

            size.scale(QtCore.QSize(rect.width()-2*margin,rect.height()-2*margin), QtCore.Qt.KeepAspectRatio)
            painter.setViewport(margin, height_doc+spacer+margin, size.width(), size.height())
            painter.setWindow(image.rect())
            painter.drawImage(0, 0, image)

            painter.restore()


            painter.save()

            penwidth = 30
            pen = QtGui.QPen(QtGui.QColor(220, 220, 220), penwidth)
            #pen = QtGui.QPen(QtGui.QColor(0, 0, 0), penwidth)
            pen.setJoinStyle(QtCore.Qt.RoundJoin);
            painter.setPen(pen)


            painter.drawRect(QtCore.QRectF(penwidth/2+5, height_doc+spacer+penwidth/2+5, size.width()+2*margin-penwidth/2-5, size.height()+2*margin-penwidth/2-5))

            #bound = QtCore.QRectF(5, height_doc+spacer+5, size.width()+2*margin-5, size.height()+2*margin-5-penwidth/2)
            #painter.drawLine(bound.bottomLeft(),bound.bottomRight())
            #painter.drawLine(bound.topLeft(),bound.topRight())

            #print(bound.bottomLeft(),bound.bottomRight())
            #print(bound.bottomLeft(),bound.bottomRight())

            plotheight = size.height()+2*margin

            painter.restore()"""

        """dialog = QPrintDialog(self.printer, self)
        if dialog.exec_():
            painter = QPainter(self.printer)
            rect = painter.viewport()
            size = self.imageLabel.pixmap().size()
            size.scale(rect.size(), Qt.KeepAspectRatio)
            painter.setViewport(rect.x(), rect.y(), size.width(), size.height())
            painter.setWindow(self.imageLabel.pixmap().rect())
            painter.drawPixmap(0, 0, self.imageLabel.pixmap())"""

        """if(dialog.exec_() != QtGui.QDialog.Accepted):
            return
        printLabel = QtGui.QLabel("Hello my printer.")
        painter = QtGui.QPainter(printer)
        printLabel.render(painter)
        painter.end()"""

    def closeEvent(self, event):
        # Kill c++ program if necessary
        self.abortCalculation()

        # Save last settings
        self.saveSettingsSystem(self.path_system_last)
        self.saveSettingsPlotter(self.path_plot_last)
        self.saveSettingsView(self.path_view_last)

        # Close everything
        super().closeEvent(event)

def main():
    app = QtGui.QApplication(sys.argv)
    form = MainWindow()
    form.show()
    rc = app.exec_()

    # Hack to avoid exit crashes
    # http://stackoverflow.com/questions/11945183/what-are-good-practices-for-avoiding-crashes-hangs-in-pyqt
    del form
    del app
    os._exit(rc)
