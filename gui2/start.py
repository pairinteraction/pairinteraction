#!/usr/bin/env python3

import sys
from pint import UnitRegistry
from pint.unit import UndefinedUnitError
from PyQt4 import QtCore, QtGui
from plotter import Ui_plotwindow # pyuic4 plotter.ui > plotter.py
import pyqtgraph as pg

import collections
from abc import ABCMeta, abstractmethod
from time import sleep
import locale
import json
import os
import multiprocessing
import subprocess
import numpy as np
from scipy import sparse
import signal
from queue import Queue
import ctypes
import sip
from scipy import constants

signal.signal(signal.SIGINT, signal.SIG_DFL)

ureg = UnitRegistry()
Q = ureg.Quantity
C = lambda s : Q(constants.value(s),constants.unit(s))

locale.setlocale(locale.LC_NUMERIC, '') # TODO auf Englisch umstellen, auch in der mit QtDesigner

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(antialias=True)

class Converter:
    converter = dict()
    dimensionless = Q('1')
    converter[str(dimensionless.dimensionality)] = dimensionless
    au_energy = C('atomic unit of energy')/C('Planck constant')
    converter[str(au_energy.dimensionality)] = au_energy
    au_bfield = C('atomic unit of mag. flux density')
    converter[str(au_bfield.dimensionality)] = au_bfield
    au_efield = C('atomic unit of electric field')
    converter[str(au_efield.dimensionality)] = au_efield
    au_length = C('atomic unit of length')
    converter[str(au_length.dimensionality)] = au_length

    @classmethod
    def toAU(self,v):
        dimensionality = v.dimensionality
        return (v/self.converter[str(dimensionality)]).to('dimensionless')

    @classmethod
    def fromAU(self, v, units):
        dimensionality = Q(1,units).dimensionality
        return (v*self.converter[str(dimensionality)]).to(units)

class Units:
    length = Q('micrometer').units
    energy = Q('gigahertz').units
    efield = Q('volt/centimeter').units
    bfield = Q('gauss').units
    angle = Q('degree').units
    dimensionless = Q('1').units    

class GUIDict(collections.MutableMapping, metaclass=ABCMeta):
    def __init__(self, ui):
        self.store = dict()
        self._setup(self.store, ui)
    
    @abstractmethod
    def _setup(self, store, ui):
        pass

    def __getitem__(self, key):
        widget = self.store[key]['widget']
        unit = self.store[key]['unit'] if 'unit' in self.store[key] else None
        
        value = None
        if isinstance(widget, QtGui.QComboBox):
            value = str(widget.currentText())
        elif isinstance(widget, QtGui.QSpinBox):
            value = int(widget.value())
        elif isinstance(widget, QtGui.QDoubleSpinBox):
            value = float(widget.value())
        elif isinstance(widget, QtGui.QLineEdit):
            try: value = locale.atof(str(widget.text()))
            except ValueError: value = None
        elif isinstance(widget, QtGui.QCheckBox):
            value = widget.checkState() == QtCore.Qt.Checked
        elif isinstance(widget, QtGui.QRadioButton):
            value = widget.isChecked()
        elif isinstance(widget, QtGui.QGroupBox):
            value = widget.isChecked()
        
        if unit is not None and value is not None:
            return Q(value, unit)
        else:
            return value

    def __setitem__(self, key, value):
        widget = self.store[key]['widget']
        unit = self.store[key]['unit'] if 'unit' in self.store[key] else None
                
        if isinstance(value,Q):
            value = value.to(unit).magnitude

        if isinstance(widget, QtGui.QComboBox):
            index = widget.findText(value)
            if index >= 0: widget.setCurrentIndex(index)
        elif isinstance(widget, QtGui.QSpinBox):
            widget.setValue(value)
        elif isinstance(widget, QtGui.QDoubleSpinBox):
            widget.setValue(value)
        elif isinstance(widget, QtGui.QLineEdit):
            if value is None: widget.setText("None")
            else: widget.setText(locale.str(value))
        elif isinstance(widget, QtGui.QCheckBox):
            if value: widget.setCheckState(QtCore.Qt.Checked)
            else: widget.setCheckState(QtCore.Qt.Unchecked)
        elif isinstance(widget, QtGui.QRadioButton):
            widget.setChecked(value)
        elif isinstance(widget, QtGui.QGroupBox):
            widget.setChecked(value)
            
    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)
    
    def load(self, f):
        params = json.load(f)
        for k, v in params.items():
            try:
                if not isinstance(v, str): raise TypeError
                v = v.replace("dimensionless","1")
                self[k] = Q(v)
            except (UndefinedUnitError, TypeError):
                self[k] = v
            
    def saveInAU(self, f, exclude = []):
        params = dict()
        for k, v in self.items():
            if k in exclude: continue
            if isinstance(v, Q): params[k] = Converter.toAU(v).magnitude                    
            else: params[k] = v
        json.dump(params, f, indent=4, sort_keys=True)
    
    def saveInOriginalunits(self, f, exclude = []):
        params = dict()
        for k, v in self.items():
            if k in exclude: continue
            if isinstance(v, Q): params[k] = str(v)
            else: params[k] = v
        json.dump(params, f, indent=4, sort_keys=True)

class SystemDict(GUIDict):
    def _setup(self, store, ui):
        store["species1"] = {'widget': ui.combobox_system_species1}
        store["species2"] = {'widget': ui.combobox_system_species2}
        store["n1"] = {'widget': ui.spinbox_system_n1, 'unit': Units.dimensionless}
        store["n2"] = {'widget': ui.spinbox_system_n2, 'unit': Units.dimensionless}
        store["l1"] = {'widget': ui.spinbox_system_l1, 'unit': Units.dimensionless}
        store["l2"] = {'widget': ui.spinbox_system_l2, 'unit': Units.dimensionless}
        store["j1"] = {'widget': ui.spinbox_system_j1, 'unit': Units.dimensionless}
        store["j2"] = {'widget': ui.spinbox_system_j2, 'unit': Units.dimensionless}
        store["m1"] = {'widget': ui.spinbox_system_m1, 'unit': Units.dimensionless}
        store["m2"] = {'widget': ui.spinbox_system_m2, 'unit': Units.dimensionless}
        store["deltaN"] = {'widget': ui.spinbox_system_deltaN, 'unit': Units.dimensionless}
        store["deltaL"] = {'widget': ui.spinbox_system_deltaL, 'unit': Units.dimensionless}
        store["deltaJ"] = {'widget': ui.spinbox_system_deltaJ, 'unit': Units.dimensionless}
        store["deltaM"] = {'widget': ui.spinbox_system_deltaM, 'unit': Units.dimensionless}
        store["deltaE1"] = {'widget': ui.lineedit_system_deltaE1, 'unit': Units.energy}
        store["deltaE2"] = {'widget': ui.lineedit_system_deltaE2, 'unit': Units.energy}
        store["deltaE"] = {'widget': ui.lineedit_system_deltaE, 'unit': Units.energy}
        store["samebasis"] = {'widget': ui.checkbox_system_samebasis}
        store["minEx"] = {'widget': ui.lineedit_system_minEx, 'unit': Units.efield}
        store["minEy"] = {'widget': ui.lineedit_system_minEy, 'unit': Units.efield}
        store["minEz"] = {'widget': ui.lineedit_system_minEz, 'unit': Units.efield}
        store["minBx"] = {'widget': ui.lineedit_system_minBx, 'unit': Units.bfield}
        store["minBy"] = {'widget': ui.lineedit_system_minBy, 'unit': Units.bfield}
        store["minBz"] = {'widget': ui.lineedit_system_minBz, 'unit': Units.bfield}
        store["maxEx"] = {'widget': ui.lineedit_system_maxEx, 'unit': Units.efield}
        store["maxEy"] = {'widget': ui.lineedit_system_maxEy, 'unit': Units.efield}
        store["maxEz"] = {'widget': ui.lineedit_system_maxEz, 'unit': Units.efield}
        store["maxBx"] = {'widget': ui.lineedit_system_maxBx, 'unit': Units.bfield}
        store["maxBy"] = {'widget': ui.lineedit_system_maxBy, 'unit': Units.bfield}
        store["maxBz"] = {'widget': ui.lineedit_system_maxBz, 'unit': Units.bfield}
        store["minR"] = {'widget': ui.lineedit_system_minR, 'unit': Units.length}
        store["maxR"] = {'widget': ui.lineedit_system_maxR, 'unit': Units.length}
        store["theta"] = {'widget': ui.lineedit_system_theta, 'unit': Units.angle}
        store["dd"] = {'widget': ui.checkbox_system_dd}
        store["dq"] = {'widget': ui.checkbox_system_dq}
        store["qq"] = {'widget': ui.checkbox_system_qq}
        store["steps"] = {'widget': ui.spinbox_system_steps, 'unit': Units.dimensionless}
        store["precision"] = {'widget': ui.lineedit_system_precision, 'unit': Units.dimensionless}
    
    def saveInAU_field1(self, f):
        self.saveInAU(f,['species2','n2','l2','m2','j2','minR','maxR','theta','dd','dq','qq']) # TODO 'deltaE2', 'deltaE'
    
    def saveInAU_field2(self, f):
        self.saveInAU(f,['species1','n1','l1','m1','j1','minR','maxR','theta','dd','dq','qq']) # TODO 'deltaE1', 'deltaE'
    
    def saveInAU_field12(self, f):
        self.saveInAU(f,['minR','maxR','theta','dd','dq','qq']) # TODO 'deltaE'
    
    def saveInAU_potential(self, f):
        self.saveInAU(f,[])

class PlotDict(GUIDict):
    def _setup(self, store, ui):
        store["minE"] = {'widget': ui.lineedit_plot_minE, 'unit': Units.energy}
        store["maxE"] = {'widget': ui.lineedit_plot_maxE, 'unit': Units.energy}
        store["lines"] = {'widget': ui.groupbox_plot_lines}
        store["points"] = {'widget': ui.groupbox_plot_points}
        store["labels"] = {'widget': ui.groupbox_plot_labels}
        store["overlap"] = {'widget': ui.groupbox_plot_overlap}
        store["szLine"] = {'widget': ui.spinbox_plot_szLine, 'unit': Units.dimensionless}
        store["szPoint"] = {'widget': ui.spinbox_plot_szPoint, 'unit': Units.dimensionless}
        store["szLabel"] = {'widget': ui.spinbox_plot_szLabel, 'unit': Units.dimensionless}
        store["transpLine"] = {'widget': ui.spinbox_plot_transpLine, 'unit': Units.dimensionless}
        store["transpPoint"] = {'widget': ui.spinbox_plot_transpPoint, 'unit': Units.dimensionless}
        store["transpLabel"] = {'widget': ui.spinbox_plot_transpLabel, 'unit': Units.dimensionless}
        store["overlapUnperturbed"] = {'widget': ui.radiobutton_plot_overlapUnperturbed}
        store["overlapDefined"] = {'widget': ui.radiobutton_plot_overlapDefined}
        store["n1"] = {'widget': ui.spinbox_plot_n1, 'unit': Units.dimensionless}
        store["n2"] = {'widget': ui.spinbox_plot_n2, 'unit': Units.dimensionless}
        store["l1"] = {'widget': ui.spinbox_plot_l1, 'unit': Units.dimensionless}
        store["l2"] = {'widget': ui.spinbox_plot_l2, 'unit': Units.dimensionless}
        store["j1"] = {'widget': ui.spinbox_plot_j1, 'unit': Units.dimensionless}
        store["j2"] = {'widget': ui.spinbox_plot_j2, 'unit': Units.dimensionless}
        store["m1"] = {'widget': ui.spinbox_plot_m1, 'unit': Units.dimensionless}
        store["m2"] = {'widget': ui.spinbox_plot_m2, 'unit': Units.dimensionless}

class Worker(QtCore.QThread):
    output = QtCore.pyqtSignal(str)

    def __init__(self, parent = None):
        super().__init__(parent)
        self.exiting = False
        self.dataqueue_field1 = Queue()
        self.dataqueue_field2 = Queue()
        self.dataqueue_field12 = Queue()
        self.dataqueue_potential = Queue()
        
    def __del__(self):
        self.exiting = True
        self.wait()
        
    def execute(self, stdout): 
        self.stdout = stdout   
        self.start()

    def run(self):
        finishedgracefully = False
        
        # Clear data queue
        with self.dataqueue_field1.mutex: self.dataqueue_field1.queue.clear()
        with self.dataqueue_field2.mutex: self.dataqueue_field2.queue.clear()
        with self.dataqueue_field12.mutex: self.dataqueue_field12.queue.clear()
        with self.dataqueue_potential.mutex: self.dataqueue_potential.queue.clear()
        
        # Parse stdout
        dim = 0
        type = 0
        current = 0
        total = 0
        
        status_type = ""
        status_progress = ""
        status_dimension = ""
        
        for line in iter(self.stdout.readline, b""):
            
            if self.exiting or not line:
                break
                
            elif line[:5] == b">>TYP":
                type = int(line[5:12].decode('utf-8'))
                status_type = ["Field map of first atom: ", "Field map of second atom: ", "Pair potential: ", "Field maps: "][type]
                status_progress = "construct matrices"
                
            elif line[:5] == b">>TOT":
                total = int(line[5:12].decode('utf-8'))
                current = 0
                
            elif line[:5] == b">>DIM":
                dim = int(line[5:12])
                status_progress = "diagonalize {} x {} matrix, {} of {} matrices already processed".format(dim, dim, current,total)
                
            elif line[:5] == b">>OUT":
                current += 1
                status_progress = "diagonalize {} x {} matrix, {} of {} matrices already processed".format(dim, dim, current,total)
                
                filenumber = int(line[5:12].decode('utf-8'))
                filename = line[13:-1].decode('utf-8')
                
                if type == 0:
                    self.dataqueue_field1.put(filename)
                elif type == 1:
                    self.dataqueue_field2.put(filename)
                elif type == 2:
                    self.dataqueue_potential.put(filename)
                elif type == 3:
                    self.dataqueue_field12.put(filename)
                    
            elif line[:5] == b">>END":
                finishedgracefully = True
                break
                
            else:
                print (line.decode('utf-8'), end="")
            
            self.output.emit(status_type + status_progress)
        
        # Clear data queue if thread has aborted
        if not finishedgracefully:
            with self.dataqueue_field1.mutex: self.dataqueue_field1.queue.clear()
            with self.dataqueue_field2.mutex: self.dataqueue_field2.queue.clear()
            with self.dataqueue_field12.mutex: self.dataqueue_field1.queue.clear()
            with self.dataqueue_potential.mutex: self.dataqueue_potential.queue.clear()
            
class BinaryLoader:
    def __init__(self):
        # types
        self.typeIds = {1008: 'int8', 1016 : 'int16', 1032 : 'int32', 1064 : 'int64', 1108 : 'uint8', 1116 : 'uint16', 1132 : 'uint32', \
            1164 : 'int64', 2032 : 'float32', 2064 : 'float64'}
        self.type_t = 'uint16'
        
        # bit masks
        self.csr_not_csc = 0x01; # xxx0: csc, xxx1: csr
        self.complex_not_real = 0x02; # xx0x: real, xx1x: complex

    def readNumber(self, f, sz = None):
        datatype = self.typeIds[np.fromfile(f, dtype=np.dtype(self.type_t), count=1)[0]]
        if sz is None: return np.fromfile(f, dtype=np.dtype(datatype), count=1)[0]
        else: return np.fromfile(f, dtype=np.dtype(datatype), count=sz)

    def readVector(self, f):
        size = self.readNumber(f)
        return self.readNumber(f, size)

    def readMatrix(self, f):
        flags = self.readNumber(f)
        rows = self.readNumber(f)
        cols = self.readNumber(f)
        if flags & self.complex_not_real: data = self.readVector(f) + self.readVector(f)*1j
        else: data = self.readVector(f)
        indices = self.readVector(f)
        indptr = np.append(self.readVector(f),len(data)) # TODO TypeError: object of type 'numpy.complex128' has no len()
        if flags & self.csr_not_csc: return sparse.csr_matrix((data, indices, indptr), shape=(rows, cols))
        else: return sparse.csc_matrix((data, indices, indptr), shape=(rows, cols))

class Eigensystem(BinaryLoader):
    def __init__(self, filename):
        super().__init__()
        
        self._filename = filename
        self._shift = 0
        
        self._params = None
        self._energies = None
        self._basis = None
        
    @property
    def params(self):
        if self._params is None:
            with open(self._filename+'.json','r') as f:
                self._params = json.load(f)
        return self._params
    
    @property
    def energies(self):
        if self._energies is None:
            with open(self._filename+'.mat','rb') as f:
                self._energies = np.real(self.readMatrix(f).diagonal())
                self._shift = f.tell()
        return self._energies
    
    @property
    def basis(self):
        if self._basis is None:
            with open(self._filename+'.mat','rb') as f:
                if self._shift > 0: f.seek(self._shift,0)
                else: self._energies = np.real(self.readMatrix(f).diagonal())
                self._basis = self.readMatrix(f)
        return self._basis

from ctypes.util import find_library
    
# see https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!msg/pyqtgraph/O-d2L6qfPoo/i1zedC2Oda4J
if sys.platform == 'win32':
    qtlib = ctypes.windll.qtgui4
    drawPoints = getattr(qtlib, '?drawPoints@QPainter@@QEAAXPEBVQPointF@@H@Z')
else:
    """
    pwd: /usr/lib/i386-linux-gnu
    list: nm -D ./libQtGui.so.4 > ~./list.txt
    """
    qtlib = ctypes.cdll.LoadLibrary(find_library("QtGui"))
    drawPoints = getattr(qtlib, '_ZN8QPainter10drawPointsEPK7QPointFi')

class PointsItem(QtGui.QGraphicsItem):
    def __init__(self, x=None, y=None, size=1, alpha=80):
        QtGui.QGraphicsItem.__init__(self)
        self.size = size
        self.alpha = alpha
        #self.pen = pg.mkPen((0,0,0,self.alpha),width=self.size,style=QtCore.Qt.CustomDashLine)
        #self.pen.setDashPattern([1, 20, 5, 4])
        self.pen = pg.mkPen((0,0,0,self.alpha),width=self.size,cosmetic=True)
        self.setData(x, y)
        #self.ItemIgnoresTransformations = True
        #self.setFlag(QtGui.QGraphicsItem.ItemIgnoresTransformations, True)
        
    def setData(self, x, y):
        if x is None:
            x = np.array([])
            y = np.array([])
        self.data = np.empty((len(x), 2), dtype=np.float)
        self.data[:,0] = x
        self.data[:,1] = y
        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()
        self.bounds = QtCore.QRectF(xmin, ymin, xmax-xmin, ymax-ymin)
        self.prepareGeometryChange()
        
        #self.qdata = [None]*len(x)
        #for n,[a,b] in enumerate(zip(x,y)):
        #    self.qdata[n] = QtCore.QPointF(a,b)
        #self.qdata = QtGui.QPolygonF(self.qdata)

    def boundingRect(self):
        return self.bounds

    def paint(self, p, *args):
        p.setPen(self.pen)
        
        #p.drawPoints(self.qdata)
        ptr = ctypes.c_void_p(sip.unwrapinstance(p))
        drawPoints(ptr, self.data.ctypes, self.data.shape[0])

# https://stackoverflow.com/questions/17103698/plotting-large-arrays-in-pyqtgraph
class MultiLine(pg.QtGui.QGraphicsPathItem):
    def __init__(self, x, y):
        """x and y are 2D arrays of shape (Nplots, Nsamples)"""
        connect = np.ones(x.shape, dtype=bool)
        connect[:,-1] = 0 # don't draw the segment between each trace
        self.path = pg.arrayToQPath(x.flatten(), y.flatten(), connect.flatten())
        pg.QtGui.QGraphicsPathItem.__init__(self, self.path)
        self.setPen(pg.mkPen((0,0,0,50),width=1))
    def shape(self): # override because QGraphicsPathItem.shape is too expensive.
        return pg.QtGui.QGraphicsItem.shape(self)
    def boundingRect(self):
        return self.path.boundingRect()

class DoublenoneValidator(QtGui.QDoubleValidator):
    def __init__(self, parent=None):
        super().__init__(parent)

    def validate(self, s, pos):
        if s == 'None':
            return (QtGui.QValidator.Acceptable, s, pos)
        
        lastpos = -1
        for c in s.lower():
            try:
                lastpos = 'none'[lastpos+1:].index(c)
            except ValueError:
                return super().validate(s, pos)
        
        return (QtGui.QValidator.Intermediate, s, pos)

    def fixup(self, s):
        return 'None'

class DoublepositiveValidator(QtGui.QDoubleValidator):
    def __init__(self, parent=None):
        super().__init__(parent)

    def validate(self, s, pos):
        status = super().validate(s, pos)
        
        if status[0] == QtGui.QValidator.Intermediate and len(s) > 0 and s[0] == '-':
            return (QtGui.QValidator.Invalid, s, pos)
        
        if status[0] == QtGui.QValidator.Acceptable and float(s) < 0:
            return (QtGui.QValidator.Invalid, s, pos)
        
        return status

    def fixup(self, s):
        return "0"

class DoubleValidator(QtGui.QDoubleValidator):
    def __init__(self, parent=None):
        super().__init__(parent)

    def validate(self, s, pos):
        return super().validate(s, pos)

    def fixup(self, s):
        return "0"

class MainWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_plotwindow()
        self.ui.setupUi(self)
        
        self.invalidQuantumnumbers = [False, False, False, False]
        
        self.samebasis_state = None
        self.samebasis = False
        
        self.systemdict = SystemDict(self.ui)
        self.plotdict = PlotDict(self.ui)
        
        self.systempath = os.getcwd()
        self.plotpath = os.getcwd()
        self.systemfile = None
        self.plotfile = None
        
        self.numprocessors = max(2,multiprocessing.cpu_count())
        self.path_base = os.path.dirname(os.path.abspath(__file__))
        self.path_workingdir = os.path.join(self.path_base,"../build/")
        self.path_cpp = os.path.join(self.path_base,"../build/pairinteraction")
        
        self.proc = None
        self.thread = Worker()
        self.timer = QtCore.QTimer()
        
        # Setup plot
        self.minE = None
        self.maxE = None
        
        self.plotOverEField = self.nonconstEField()
        
        for plotarea in [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot, self.ui.graphicsview_potential_plot]:
            plotarea.setDownsampling(ds=True, auto=True, mode='peak')
            plotarea.setClipToView(True)
            plotarea.setLabel('left', 'Energy ('+str(Units.energy)+')')
        
        if self.plotOverEField:
            for plotarea in [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot]: plotarea.setLabel('bottom', 'Electric field strength ('+str(Units.efield)+')')
        else:
            for plotarea in [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot]: plotarea.setLabel('bottom', 'Magnetic field strength ('+str(Units.bfield)+')')
        self.ui.graphicsview_potential_plot.setLabel('bottom', 'Interatomic distance ('+str(Units.length)+')')
        
        # Set validators
        validator_double = DoubleValidator()
        validator_doublenone = DoublenoneValidator()
        validator_doublepositive = DoublepositiveValidator()
        self.ui.lineedit_system_deltaE1.setValidator(validator_doublepositive)
        self.ui.lineedit_system_deltaE2.setValidator(validator_doublepositive)
        self.ui.lineedit_system_deltaE.setValidator(validator_doublepositive)
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
        self.ui.lineedit_system_precision.setValidator(validator_doublepositive)
        self.ui.lineedit_plot_minE.setValidator(validator_doublenone)
        self.ui.lineedit_plot_maxE.setValidator(validator_doublenone)
        
        # Connect signals and slots
        """self.ui.lineedit_system_deltaE1.textChanged.connect(self.checkState)
        self.ui.lineedit_system_deltaE2.textChanged.connect(self.checkState)
        self.ui.lineedit_system_deltaE.textChanged.connect(self.checkState)
        self.ui.lineedit_system_minEx.textChanged.connect(self.checkState)
        self.ui.lineedit_system_minEy.textChanged.connect(self.checkState)
        self.ui.lineedit_system_minEz.textChanged.connect(self.checkState)
        self.ui.lineedit_system_maxEx.textChanged.connect(self.checkState)
        self.ui.lineedit_system_maxEy.textChanged.connect(self.checkState)
        self.ui.lineedit_system_maxEz.textChanged.connect(self.checkState)
        self.ui.lineedit_system_minBx.textChanged.connect(self.checkState)
        self.ui.lineedit_system_minBy.textChanged.connect(self.checkState)
        self.ui.lineedit_system_minBz.textChanged.connect(self.checkState)
        self.ui.lineedit_system_maxBx.textChanged.connect(self.checkState)
        self.ui.lineedit_system_maxBy.textChanged.connect(self.checkState)
        self.ui.lineedit_system_maxBz.textChanged.connect(self.checkState)
        self.ui.lineedit_system_minR.textChanged.connect(self.checkState)
        self.ui.lineedit_system_minR.textChanged.connect(self.checkState)
        self.ui.lineedit_system_maxR.textChanged.connect(self.checkState)
        self.ui.lineedit_system_theta.textChanged.connect(self.checkState)
        self.ui.lineedit_system_precision.textChanged.connect(self.checkState)
        self.ui.lineedit_plot_minE.textChanged.connect(self.checkState)
        self.ui.lineedit_plot_maxE.textChanged.connect(self.checkState)"""
        
        self.ui.spinbox_system_n1.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_system_n2.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_system_l1.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_system_l2.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_system_j1.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_system_j2.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_system_m1.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_system_m2.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_plot_n1.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_plot_n2.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_plot_l1.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_plot_l2.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_plot_j1.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_plot_j2.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_plot_m1.valueChanged.connect(self.validateQuantumnumbers)
        self.ui.spinbox_plot_m2.valueChanged.connect(self.validateQuantumnumbers)

        self.ui.spinbox_system_j1.editingFinished.connect(self.validateHalfinteger)
        self.ui.spinbox_system_j2.editingFinished.connect(self.validateHalfinteger)
        self.ui.spinbox_system_m1.editingFinished.connect(self.validateHalfinteger)
        self.ui.spinbox_system_m2.editingFinished.connect(self.validateHalfinteger)
        self.ui.spinbox_plot_j1.editingFinished.connect(self.validateHalfinteger)
        self.ui.spinbox_plot_j2.editingFinished.connect(self.validateHalfinteger)
        self.ui.spinbox_plot_m1.editingFinished.connect(self.validateHalfinteger)
        self.ui.spinbox_plot_m2.editingFinished.connect(self.validateHalfinteger)
        
        self.ui.lineedit_system_deltaE1.textChanged.connect(self.forbidSamebasis)
        self.ui.lineedit_system_deltaE2.textChanged.connect(self.forbidSamebasis)
        self.ui.combobox_system_species1.currentIndexChanged.connect(self.forbidSamebasis)
        self.ui.combobox_system_species2.currentIndexChanged.connect(self.forbidSamebasis)
        
        self.ui.action_system_open.triggered.connect(self.openSystemConf)
        self.ui.action_system_save.triggered.connect(self.saveSystemConf)
        self.ui.action_plot_open.triggered.connect(self.openPlotConf)
        self.ui.action_plot_save.triggered.connect(self.savePlotConf)
        self.ui.pushbutton_field1_calc.clicked.connect(self.startCalc)
        self.ui.pushbutton_field2_calc.clicked.connect(self.startCalc)
        self.ui.pushbutton_potential_calc.clicked.connect(self.startCalc)
        self.thread.output.connect(self.updateStatus)
        self.timer.timeout.connect(self.checkForData)
        
        # Load last settings
        self.path_system_last = os.path.join(self.path_base,"lastsystem.json")
        self.path_plot_last = os.path.join(self.path_base,"lastplotter.json")
        
        if os.path.isfile(self.path_system_last):
            with open(self.path_system_last, 'r') as f:
                self.systemdict.load(f)
        
        if os.path.isfile(self.path_plot_last):
            with open(self.path_plot_last, 'r') as f:
                self.plotdict.load(f)
        
        # Emit change-signals in order to let the validation run
        self.ui.spinbox_system_n1.valueChanged.emit(self.ui.spinbox_system_n1.value())
        self.ui.spinbox_system_n2.valueChanged.emit(self.ui.spinbox_system_n2.value())
        self.ui.spinbox_system_l1.valueChanged.emit(self.ui.spinbox_system_l1.value())
        self.ui.spinbox_system_l2.valueChanged.emit(self.ui.spinbox_system_l2.value())
        self.ui.spinbox_system_j1.valueChanged.emit(self.ui.spinbox_system_j1.value())
        self.ui.spinbox_system_j2.valueChanged.emit(self.ui.spinbox_system_j2.value())
        self.ui.spinbox_system_m1.valueChanged.emit(self.ui.spinbox_system_m1.value())
        self.ui.spinbox_system_m2.valueChanged.emit(self.ui.spinbox_system_m2.value())
        self.ui.spinbox_plot_n1.valueChanged.emit(self.ui.spinbox_plot_n1.value())
        self.ui.spinbox_plot_n2.valueChanged.emit(self.ui.spinbox_plot_n2.value())
        self.ui.spinbox_plot_l1.valueChanged.emit(self.ui.spinbox_plot_l1.value())
        self.ui.spinbox_plot_l2.valueChanged.emit(self.ui.spinbox_plot_l2.value())
        self.ui.spinbox_plot_j1.valueChanged.emit(self.ui.spinbox_plot_j1.value())
        self.ui.spinbox_plot_j2.valueChanged.emit(self.ui.spinbox_plot_j2.value())
        self.ui.spinbox_plot_m1.valueChanged.emit(self.ui.spinbox_plot_m1.value())
        self.ui.spinbox_plot_m2.valueChanged.emit(self.ui.spinbox_plot_m2.value())
        self.ui.lineedit_system_deltaE1.textChanged.emit(self.ui.lineedit_system_deltaE1.text())
        self.ui.lineedit_system_deltaE2.textChanged.emit(self.ui.lineedit_system_deltaE2.text())
        self.ui.combobox_system_species1.currentIndexChanged.emit(self.ui.combobox_system_species1.currentIndex())
        self.ui.combobox_system_species2.currentIndexChanged.emit(self.ui.combobox_system_species2.currentIndex())
    
    def resizeEvent(self, event):
        super().resizeEvent(event)
        if sys.platform == "darwin": QtGui.QApplication.processEvents() # hack to circumvent the no-redraw-after-resizing-bug
    
    def nonconstEField(self):
        minEField = np.linalg.norm([self.systemdict['minEx'].magnitude,self.systemdict['minEy'].magnitude,self.systemdict['minEz'].magnitude]) # TODO
        maxEField = np.linalg.norm([self.systemdict['maxEx'].magnitude,self.systemdict['maxEy'].magnitude,self.systemdict['maxEz'].magnitude]) # TODO
        minBField = np.linalg.norm([self.systemdict['minBx'].magnitude,self.systemdict['minBy'].magnitude,self.systemdict['minBz'].magnitude]) # TODO
        maxBField = np.linalg.norm([self.systemdict['maxBx'].magnitude,self.systemdict['maxBy'].magnitude,self.systemdict['maxBz'].magnitude]) # TODO
        return (maxEField != minEField) or (maxBField == minBField)        
    
    def checkForData(self):

        # check if there is some new data to plot
        if not self.thread.dataqueue_potential.empty():
        
            # Draw hole buffer at once, at least if it is no very large
            x = np.array([])
            y = np.array([])
            
            while not self.thread.dataqueue_potential.empty():
                eigensystem = Eigensystem(self.thread.dataqueue_potential.get())
                energies = eigensystem.energies
                position = float(eigensystem.params["R"])
                
                if self.minE is not None: energies = energies[energies >= self.minE]
                if self.maxE is not None: energies = energies[energies <= self.maxE]
                                
                x = np.append(x,position*np.ones_like(energies))
                y = np.append(y,energies)
                
                if len(x) > 5000: break                
            
            if len(x) > 0:
                x *= Converter.fromAU(1,Units.length).magnitude
                y *= Converter.fromAU(1,Units.energy).magnitude
            
                if self.ui.groupbox_plot_points.isChecked():
                    size = self.ui.spinbox_plot_szPoint.value()
                    alpha = self.ui.spinbox_plot_transpPoint.value()*255
                    curve = PointsItem(x, y, size, alpha)
                    self.ui.graphicsview_potential_plot.addItem(curve)
                
                self.ui.graphicsview_potential_plot.repaint()
        
        # check if there is some new data to plot
        if not self.thread.dataqueue_field1.empty():
        
            # Draw hole buffer at once, at least if it is no very large
            x = np.array([])
            y = np.array([])
            
            while not self.thread.dataqueue_field1.empty():
                eigensystem = Eigensystem(self.thread.dataqueue_field1.get())
                energies = eigensystem.energies
                if self.plotOverEField: position = np.linalg.norm([float(eigensystem.params["Ex"]),float(eigensystem.params["Ey"]),float(eigensystem.params["Ez"])])  # TODO
                else: position = np.linalg.norm([float(eigensystem.params["Bx"]),float(eigensystem.params["By"]),float(eigensystem.params["Bz"])])  # TODO
                
                if self.minE is not None: energies = energies[energies >= self.minE]
                if self.maxE is not None: energies = energies[energies <= self.maxE]
                
                x = np.append(x,position*np.ones_like(energies))
                y = np.append(y,energies)
                
                if len(x) > 5000: break                
            
            if len(x) > 0:
                if self.plotOverEField: x *= Converter.fromAU(1,Units.efield).magnitude
                else: x *= Converter.fromAU(1,Units.bfield).magnitude
                y *= Converter.fromAU(1,Units.energy).magnitude
            
                if self.ui.groupbox_plot_points.isChecked():
                    size = self.ui.spinbox_plot_szPoint.value()
                    alpha = self.ui.spinbox_plot_transpPoint.value()*255
                    curve = PointsItem(x, y, size, alpha)
                    self.ui.graphicsview_field1_plot.addItem(curve)
                
                self.ui.graphicsview_field1_plot.repaint()
        
        # check if there is some new data to plot
        if not self.thread.dataqueue_field2.empty():
        
            # Draw hole buffer at once, at least if it is no very large
            x = np.array([])
            y = np.array([])
            
            while not self.thread.dataqueue_field2.empty():
                eigensystem = Eigensystem(self.thread.dataqueue_field2.get())
                energies = eigensystem.energies
                if self.plotOverEField: position = np.linalg.norm([float(eigensystem.params["Ex"]),float(eigensystem.params["Ey"]),float(eigensystem.params["Ez"])])  # TODO
                else: position = np.linalg.norm([float(eigensystem.params["Bx"]),float(eigensystem.params["By"]),float(eigensystem.params["Bz"])])  # TODO
                
                if self.minE is not None: energies = energies[energies >= self.minE]
                if self.maxE is not None: energies = energies[energies <= self.maxE]
                
                x = np.append(x,position*np.ones_like(energies))
                y = np.append(y,energies)
                
                if len(x) > 5000: break                
            
            if len(x) > 0:
                if self.plotOverEField: x *= Converter.fromAU(1,Units.efield).magnitude
                else: x *= Converter.fromAU(1,Units.bfield).magnitude
                y *= Converter.fromAU(1,Units.energy).magnitude
            
                if self.ui.groupbox_plot_points.isChecked():
                    size = self.ui.spinbox_plot_szPoint.value()
                    alpha = self.ui.spinbox_plot_transpPoint.value()*255
                    curve = PointsItem(x, y, size, alpha)
                    self.ui.graphicsview_field2_plot.addItem(curve)
                
                self.ui.graphicsview_field2_plot.repaint()
        
        # check if there is some new data to plot
        if not self.thread.dataqueue_field12.empty():
        
            # Draw hole buffer at once, at least if it is no very large
            x = np.array([])
            y = np.array([])
            
            while not self.thread.dataqueue_field12.empty():
                eigensystem = Eigensystem(self.thread.dataqueue_field12.get())
                energies = eigensystem.energies
                if self.plotOverEField: position = np.linalg.norm([float(eigensystem.params["Ex"]),float(eigensystem.params["Ey"]),float(eigensystem.params["Ez"])])  # TODO
                else: position = np.linalg.norm([float(eigensystem.params["Bx"]),float(eigensystem.params["By"]),float(eigensystem.params["Bz"])])  # TODO
                
                if self.minE is not None: energies = energies[energies >= self.minE]
                if self.maxE is not None: energies = energies[energies <= self.maxE]
                
                x = np.append(x,position*np.ones_like(energies))
                y = np.append(y,energies)
                
                if len(x) > 5000: break                
            
            if len(x) > 0:
                if self.plotOverEField: x *= Converter.fromAU(1,Units.efield).magnitude
                else: x *= Converter.fromAU(1,Units.bfield).magnitude
                y *= Converter.fromAU(1,Units.energy).magnitude
            
                if self.ui.groupbox_plot_points.isChecked():
                    size = self.ui.spinbox_plot_szPoint.value()
                    alpha = self.ui.spinbox_plot_transpPoint.value()*255
                    curve1 = PointsItem(x, y, size, alpha)
                    curve2 = PointsItem(x, y, size, alpha)
                    self.ui.graphicsview_field1_plot.addItem(curve1)
                    self.ui.graphicsview_field2_plot.addItem(curve2)
                
                self.ui.graphicsview_field1_plot.repaint()
                self.ui.graphicsview_field2_plot.repaint()
            
            
            
            
            
            #Converter.fromAU(1,'gauss').magnitude
            #Converter.fromAU(1,'volt/centimeter').magnitude
            
            
        
            #self.ui.graphicsview_potential_plot.clear()
        
            #float(eigensystem.params["Ez"])*np.ones_like(eigensystem.energies)
            """curve = pg.ScatterPlotItem(x=x, y=y, 
                               pen='w', brush='b', size=3, 
                               pxMode=True,antialias=True)"""
            """curve = gl.GLScatterPlotItem(pos=np.transpose([x,y,np.zeros_like(x)]))"""
            
            
            
            # TODO lines
            #self.ui.graphicsview_potential_plot.plot(x,y)
        
            """self.ui.graphicsview_potential_plot.plot(float(eigensystem.params["R"])*np.ones_like(eigensystem.energies), eigensystem.energies, pen=None, symbol='o')
            """
            #self.ui.graphicsview_potential_plot.repaint()
        
            #self.ui.graphicsview_potential_plot.plot(eigensystem.energies)
        
        # check if thread has finished
        if self.thread.isFinished() and self.thread.dataqueue_field1.empty() and self.thread.dataqueue_field2.empty() and self.thread.dataqueue_field12.empty() and self.thread.dataqueue_potential.empty():
        
            # Enable auto range again
            #self.ui.graphicsview_potential_plot.autoRange()
        
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
            self.senderbutton.setText("Start calculation")
            
            # Reset status bar
            self.ui.statusbar.showMessage('')
    
    # https://groups.google.com/forum/#!msg/pyqtgraph/srQqVW9bqPg/CuCgyxzWo14J
    # TODO
    
    """@QtCore.pyqtSlot()
    def checkState(self):
        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = '#ffffff'
        elif state == QtGui.QValidator.Intermediate:
            color = '#fff79a'
        else:
            color = '#f6989d'
        sender.setStyleSheet('QLineEdit { background-color: %s }' % color)"""
    
    @QtCore.pyqtSlot(str)
    def forbidSamebasis(self):
        if self.ui.lineedit_system_deltaE1.text() != self.ui.lineedit_system_deltaE2.text() or \
            self.ui.combobox_system_species1.currentIndex() != self.ui.combobox_system_species2.currentIndex():
            self.ui.checkbox_system_samebasis.setEnabled(False)
            if self.samebasis_state is None:
                self.samebasis_state = self.ui.checkbox_system_samebasis.checkState()
                self.ui.checkbox_system_samebasis.setCheckState(QtCore.Qt.Unchecked)
        else:
            self.ui.checkbox_system_samebasis.setEnabled(True)
            if self.samebasis_state is not None:
                self.ui.checkbox_system_samebasis.setCheckState(self.samebasis_state)
                self.samebasis_state = None
    
    @QtCore.pyqtSlot(str)   
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
            
                if n.value()-1 < l.value():
                    n_err |= True
                    l_err |= True
            
                if l.value()+0.5 < j.value():
                    l_err |= True
                    j_err |= True
            
                if j.value() < abs(m.value()):
                    j_err |= True
                    m_err |= True
            
                """if n_err: n.setStyleSheet('QSpinBox { background-color: #f6989d }')
                else: n.setStyleSheet('QSpinBox { background-color: #ffffff }')
            
                if l_err: l.setStyleSheet('QSpinBox { background-color: #f6989d }')
                else: l.setStyleSheet('QSpinBox { background-color: #ffffff }')
            
                if j_err: j.setStyleSheet('QDoubleSpinBox { background-color: #f6989d }')
                else: j.setStyleSheet('QDoubleSpinBox { background-color: #ffffff }')
            
                if m_err: m.setStyleSheet('QDoubleSpinBox { background-color: #f6989d }')
                else: m.setStyleSheet('QDoubleSpinBox { background-color: #ffffff }')"""
                
                if n_err or l_err or j_err or m_err:
                    self.invalidQuantumnumbers[i] = True
                else:
                    self.invalidQuantumnumbers[i] = False
        
        if np.any(self.invalidQuantumnumbers):
            self.ui.statusbar.showMessage('Invalide quantum numbers specified.')
            """self.ui.pushbutton_field1_calc.setEnabled(False)
            self.ui.pushbutton_field2_calc.setEnabled(False)
            self.ui.pushbutton_potential_calc.setEnabled(False)"""
        else:
            self.ui.statusbar.showMessage('')
            """self.ui.pushbutton_field1_calc.setEnabled(True)
            self.ui.pushbutton_field2_calc.setEnabled(True)
            self.ui.pushbutton_potential_calc.setEnabled(True)"""
    
    @QtCore.pyqtSlot(str)   
    def validateHalfinteger(self):
        value = self.sender().value()
        self.sender().setValue(np.floor(value)+0.5)
    
    @QtCore.pyqtSlot(str)      
    def updateStatus(self, msg):
        if not self.proc is None: self.ui.statusbar.showMessage(msg)
            
    @QtCore.pyqtSlot()
    def startCalc(self):   
        if self.proc is None:
            if np.any(self.invalidQuantumnumbers):
                QtGui.QMessageBox.critical(self, "Message", "Invalide quantum numbers specified.")
                
            else:
                path_config = os.path.join(self.path_base,"../build/system.json")
        
                # Change buttons
                self.senderbutton = self.sender()
                if self.senderbutton != self.ui.pushbutton_field1_calc:
                    self.ui.pushbutton_field1_calc.setEnabled(False)
                if self.senderbutton != self.ui.pushbutton_field2_calc:
                    self.ui.pushbutton_field2_calc.setEnabled(False)
                if self.senderbutton != self.ui.pushbutton_potential_calc:
                    self.ui.pushbutton_potential_calc.setEnabled(False)
                self.senderbutton.setText("Abort calculation")
                
                # Store, whether the same basis should be used for both atoms
                self.samebasis = self.ui.checkbox_system_samebasis.checkState() == QtCore.Qt.Checked
            
                # Clear plots and set them up # disable auto range (for higher update speed) # TODO set axis range manually, self.ui.graphicsview_field1_plot.disableAutoRange()
                self.plotOverEField = self.nonconstEField()
            
                if self.senderbutton in [self.ui.pushbutton_field1_calc, self.ui.pushbutton_potential_calc] or self.samebasis:
                    self.ui.graphicsview_field1_plot.clear()
                    self.ui.graphicsview_field1_plot.setLabel('left', 'Energy ('+str(Units.energy)+')')
                    if self.plotOverEField: self.ui.graphicsview_field1_plot.setLabel('bottom', 'Electric field ('+str(Units.efield)+')')
                    else: self.ui.graphicsview_field1_plot.setLabel('bottom', 'Magnetic field ('+str(Units.bfield)+')')
                if self.senderbutton in [self.ui.pushbutton_field2_calc, self.ui.pushbutton_potential_calc] or self.samebasis:
                    self.ui.graphicsview_field2_plot.clear()
                    self.ui.graphicsview_field2_plot.setLabel('left', 'Energy ('+str(Units.energy)+')')
                    if self.plotOverEField: self.ui.graphicsview_field2_plot.setLabel('bottom', 'Electric field ('+str(Units.efield)+')')
                    else: self.ui.graphicsview_field2_plot.setLabel('bottom', 'Magnetic field ('+str(Units.bfield)+')')
                if self.senderbutton == self.ui.pushbutton_potential_calc:
                    self.ui.graphicsview_potential_plot.clear()
                    self.ui.graphicsview_potential_plot.setLabel('left', 'Energy ('+str(Units.energy)+')')
                    self.ui.graphicsview_potential_plot.setLabel('bottom', 'Interatomic distance ('+str(Units.length)+')')
            
                # Set limits
                self.minE = self.plotdict["minE"]
                if self.minE is not None:
                    self.minE = Converter.toAU(self.minE).magnitude
                
                self.maxE = self.plotdict["maxE"]
                if self.maxE is not None:
                    self.maxE = Converter.toAU(self.maxE).magnitude
            
                # Save configuration to json file
                with open(path_config, 'w') as f:
                    if self.senderbutton == self.ui.pushbutton_field1_calc:
                        if self.samebasis:
                            self.systemdict.saveInAU_field12(f)
                        else:
                            self.systemdict.saveInAU_field1(f)
                    elif self.senderbutton == self.ui.pushbutton_field2_calc:
                        if self.samebasis:
                            self.systemdict.saveInAU_field12(f)
                        else:
                            self.systemdict.saveInAU_field2(f)
                    elif self.senderbutton == self.ui.pushbutton_potential_calc:
                        self.systemdict.saveInAU_potential(f)
                    
                # Start c++ process # TODO use pairinteraction-real if possible
                self.proc = subprocess.Popen(["mpiexec","-n","%d"%self.numprocessors,self.path_cpp,"-c",path_config],
                    stdout=subprocess.PIPE, cwd=self.path_workingdir)
        
                # Start thread that collects the output
                self.thread.execute(self.proc.stdout)
            
                # Start timer used for processing the results
                self.timer.start(0)
            
        else:
            # Kill c++ process - this terminates the self.thread and clears self.thread.dataqueue, too
            self.proc.kill()            
    
    @QtCore.pyqtSlot()
    def saveSystemConf(self):
        path = self.systemfile if self.systemfile is not None else self.systempath
        filename = QtGui.QFileDialog.getSaveFileName(self, \
            "Save system configuration",path, "json (*.json)")
        
        if filename:
            with open(filename, 'w') as f:
                self.systemdict.saveInOriginalunits(f)
            self.systemfile = filename
            self.systempath = os.path.dirname(filename)
    
    @QtCore.pyqtSlot()
    def savePlotConf(self):
        path = self.plotfile if self.plotfile is not None else self.plotpath
        filename = QtGui.QFileDialog.getSaveFileName(self, \
            "Save plot configuration",path, "json (*.json)")
        
        if filename:
            with open(filename, 'w') as f:
                self.plotdict.saveInOriginalunits(f)
            self.plotfile = filename
            self.plotpath = os.path.dirname(filename)
    
    @QtCore.pyqtSlot()
    def openSystemConf(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, \
            "Open system configuration",self.systempath, "json (*.json)")
        
        if filename:
            with open(filename, 'r') as f:
                self.systemdict.load(f)
            self.systemfile = filename
            self.systempath = os.path.dirname(filename)
    
    @QtCore.pyqtSlot()
    def openPlotConf(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, \
            "Open plot configuration",self.systempath, "json (*.json)")
        
        if not (filename == ""):
            with open(filename, 'r') as f:
                self.plotdict.load(f)
            self.plotfile = filename
            self.plotpath = os.path.dirname(filename)
            
    def closeEvent(self, event):
        # Kill c++ program if necessary 
        if self.proc is not None: self.proc.kill()
        
        # Save last settings
        self.path_system_last = os.path.join(self.path_base,"lastsystem.json")
        self.path_plot_last = os.path.join(self.path_base,"lastplotter.json")
        
        with open(self.path_system_last, 'w') as f:
            self.systemdict.saveInOriginalunits(f)
        
        with open(self.path_plot_last, 'w') as f:
            self.plotdict.saveInOriginalunits(f)
        
        # Close everything
        super().closeEvent(event)

def main():
    app = QtGui.QApplication(sys.argv)
    form = MainWindow()
    form.show()
    app.exec_()
    
if __name__ == "__main__":
    main()
