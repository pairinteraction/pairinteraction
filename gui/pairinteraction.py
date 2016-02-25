#!/usr/bin/env python3

import sys
from pint import UnitRegistry
from pint.unit import UndefinedUnitError
from PyQt4 import QtCore, QtGui
from plotter import Ui_plotwindow # pyuic4 plotter.ui > plotter.py or py3uic4 plotter.ui > plotter.py
import pyqtgraph as pg

import collections
from abc import ABCMeta, abstractmethod
from time import sleep, time
from datetime import timedelta
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
import psutil

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
        
    def __getattr__(self, key):
        return self.__getitem__(key)

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
            
    def paramsInAU(self, f, exclude = []):
        params = dict()
        for k, v in self.items():
            if k in exclude: continue
            if isinstance(v, Q): params[k] = Converter.toAU(v).magnitude                    
            else: params[k] = v
        return params
    
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
        store["deltaNSingle"] = {'widget': ui.spinbox_system_deltaNSingle, 'unit': Units.dimensionless}
        store["deltaLSingle"] = {'widget': ui.spinbox_system_deltaLSingle, 'unit': Units.dimensionless}
        store["deltaJSingle"] = {'widget': ui.spinbox_system_deltaJSingle, 'unit': Units.dimensionless}
        store["deltaMSingle"] = {'widget': ui.spinbox_system_deltaMSingle, 'unit': Units.dimensionless}
        store["deltaNPair"] = {'widget': ui.spinbox_system_deltaNPair, 'unit': Units.dimensionless}
        store["deltaLPair"] = {'widget': ui.spinbox_system_deltaLPair, 'unit': Units.dimensionless}
        store["deltaJPair"] = {'widget': ui.spinbox_system_deltaJPair, 'unit': Units.dimensionless}
        store["deltaMPair"] = {'widget': ui.spinbox_system_deltaMPair, 'unit': Units.dimensionless}
        store["deltaESingle"] = {'widget': ui.lineedit_system_deltaESingle, 'unit': Units.energy}
        store["deltaEPair"] = {'widget': ui.lineedit_system_deltaEPair, 'unit': Units.energy}
        store["pairbasisSame"] = {'widget': ui.radiobutton_system_pairbasisSame}
        store["pairbasisDefined"] = {'widget': ui.radiobutton_system_pairbasisDefined}
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
        params = self.paramsInAU(f,['pairbasisSame','pairbasisDefined','species2','n2','l2','m2','j2','minR','maxR','theta','dd','dq','qq','deltaNPair','deltaLPair','deltaMPair','deltaJPair','deltaEPair'])
        json.dump(params, f, indent=4, sort_keys=True)
    
    def saveInAU_field2(self, f):
        params = self.paramsInAU(f,['pairbasisSame','pairbasisDefined','species1','n1','l1','m1','j1','minR','maxR','theta','dd','dq','qq','deltaNPair','deltaLPair','deltaMPair','deltaJPair','deltaEPair'])
        json.dump(params, f, indent=4, sort_keys=True)
    
    def saveInAU_field12(self, f):
        params = self.paramsInAU(f,['pairbasisSame','pairbasisDefined','minR','maxR','theta','dd','dq','qq','deltaNPair','deltaLPair','deltaMPair','deltaJPair','deltaEPair'])
        json.dump(params, f, indent=4, sort_keys=True)
    
    def saveInAU_potential(self, f):
        params = self.paramsInAU(f,['pairbasisSame','pairbasisDefined'])
        if self["pairbasisSame"]:
            params["deltaNPair"] = -1
            params["deltaLPair"] = -1
            params["deltaJPair"] = -1
            params["deltaMPair"] = -1
        json.dump(params, f, indent=4, sort_keys=True)

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

    def __init__(self, parent = None):
        super().__init__(parent)
        self.exiting = False
        self.samebasis = False
        self.message = ""
        self.basisfile_field1 = ""
        self.basisfile_field2 = ""
        self.basisfile_potential = ""
        self.dataqueue_field1 = Queue()
        self.dataqueue_field2 = Queue()
        self.dataqueue_potential = Queue()
        
    def __del__(self):
        self.exiting = True
        self.wait()
        
    def execute(self, stdout): 
        self.stdout = stdout   
        self.start()
        
    def clear(self):
        with self.dataqueue_field1.mutex: self.dataqueue_field1.queue.clear()
        with self.dataqueue_field2.mutex: self.dataqueue_field2.queue.clear()
        with self.dataqueue_potential.mutex: self.dataqueue_potential.queue.clear()

    def run(self):
        finishedgracefully = False
                
        self.message = ""
        
        # Clear filenames
        self.basisfile_field1 = ""
        self.basisfile_field2 = ""
        self.basisfile_potential = ""
        
        # Clear data queue
        with self.dataqueue_field1.mutex: self.dataqueue_field1.queue.clear()
        with self.dataqueue_field2.mutex: self.dataqueue_field2.queue.clear()
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
                
                if type == 3: self.samebasis = True
                elif type == 0 or type == 1: self.samebasis = False
            
            elif line[:5] == b">>BAS":
                basissize = int(line[5:12].decode('utf-8'))
                status_progress = "construct matrices using {} basis vectors".format(basissize)
            
            elif line[:5] == b">>STA":
                filename = line[6:-1].decode('utf-8')
                if type == 0 or type == 3:
                    self.basisfile_field1 = filename
                elif type == 1:
                    self.basisfile_field2 = filename
                elif type == 2:
                    self.basisfile_potential = filename
                
            elif line[:5] == b">>TOT":
                total = int(line[5:12].decode('utf-8'))
                current = 0
                
            elif line[:5] == b">>DIM":
                dim = int(line[5:12])
                status_progress = "diagonalize {} x {} matrix, {} of {} matrices processed".format(dim, dim, current,total)
                
            elif line[:5] == b">>OUT":
                current += 1
                status_progress = "diagonalize {} x {} matrix, {} of {} matrices processed".format(dim, dim, current,total)
                
                filenumber = int(line[5:12].decode('utf-8'))
                filestep = int(line[12:19].decode('utf-8'))
                filename = line[20:-1].decode('utf-8')
                
                if type == 0 or type == 3:
                    self.dataqueue_field1.put([filestep,filename])
                elif type == 1:
                    self.dataqueue_field2.put([filestep,filename])
                elif type == 2:
                    self.dataqueue_potential.put([filestep,filename])
                    
            elif line[:5] == b">>END":
                finishedgracefully = True
                break
                
            else:
                print (line.decode('utf-8'), end="")
            
            self.message = status_type + status_progress
        
        # Clear data queue if thread has aborted
        if not finishedgracefully:
            self.clear()
            
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
        indptr = np.append(self.readVector(f),len(data))
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
    def __init__(self, x=None, y=None, size=1, alpha=80, color=(0,0,0)):
        QtGui.QGraphicsItem.__init__(self)
        self.size = size
        self.alpha = alpha
        #self.pen = pg.mkPen((0,0,0,self.alpha),width=self.size,style=QtCore.Qt.CustomDashLine)
        #self.pen.setDashPattern([1, 20, 5, 4])
        self.pen = pg.mkPen(color+(self.alpha,),width=self.size,cosmetic=True)
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
        
        if status[0] == QtGui.QValidator.Acceptable and locale.atof(s) < 0:
            return (QtGui.QValidator.Invalid, s, pos)
        
        return status

    def fixup(self, s):
        return "0"

class DoubledeltaValidator(QtGui.QDoubleValidator):
    def __init__(self, parent=None):
        super().__init__(parent)

    def validate(self, s, pos):
        status = super().validate(s, pos)
        
        if status[0] == QtGui.QValidator.Acceptable and locale.atof(s) < 0 and locale.atof(s) != -1:
            return (QtGui.QValidator.Intermediate, s, pos)
        
        return status

    def fixup(self, s):
        if locale.atof(s) < 0: return "-1"
        return "0"

class DoubleValidator(QtGui.QDoubleValidator):
    def __init__(self, parent=None):
        super().__init__(parent)

    def validate(self, s, pos):
        return super().validate(s, pos)

    def fixup(self, s):
        return "0"

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

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
        self.path_base = os.path.dirname(os.path.realpath(__file__))
        self.path_workingdir = os.path.join(self.path_base,"../calc/")
        self.path_cpp_real = os.path.join(self.path_base,"../calc/pairinteraction-real")
        self.path_cpp_complex = os.path.join(self.path_base,"../calc/pairinteraction-complex")
        if os.name == 'nt': self.path_out = os.path.join(os.path.expanduser('~user'), "pairinteraction/")
        else: self.path_out = os.path.join(os.path.expanduser('~'), ".pairinteraction/")
        self.path_system_last = os.path.join(self.path_out,"lastsystem.json")
        self.path_plot_last = os.path.join(self.path_out,"lastplotter.json")
        self.path_config = os.path.join(self.path_out,"conf.json")
        
        self.proc = None
        
        self.thread = Worker()
        self.timer = QtCore.QTimer()
        
        self.momentumcolors = [(55,126,184),(77,175,74),(228,26,28),(152,78,163),(0,0,0),(255//5,255//5,255//5)] # s, p, d, f, other, undetermined
        
        self.momentummat = [None]*3
        self.labelmat = [None]*3
        self.labelstates = [None]*3
        self.momentumstrings = [None]*3
        
        # TODOs
        self.ui.groupbox_plot_lines.setEnabled(False)
        self.ui.groupbox_plot_overlap.setEnabled(False)
        self.ui.lineedit_system_theta.setEnabled(False)
        self.ui.lineedit_system_precision.setEnabled(False)
        self.ui.pushbutton_field1_save.setEnabled(False)
        self.ui.pushbutton_field2_save.setEnabled(False)
        self.ui.pushbutton_potential_save.setEnabled(False)
        self.ui.checkbox_system_dq.setEnabled(False)
        self.ui.checkbox_system_qq.setEnabled(False)
        self.ui.lineedit_system_deltaESingle.setStatusTip('A value of -1 means that there are no restrictions for single atom energies.')
        self.ui.lineedit_system_deltaEPair.setStatusTip('A value of -1 means that there are no restrictions for pair energies.')
        self.ui.radiobutton_system_pairbasisDefined.setText("use the restrictions below")
        
        # Create directories
        if not os.path.exists(self.path_out):
            os.makedirs(self.path_out)	
            if os.name == 'nt':
                ret = ctypes.windll.kernel32.SetFileAttributesW(self.path_out,FILE_ATTRIBUTE_HIDDEN)
                if not ret: raise ctypes.WinError()
                
        # Load last settings
        try: # TODO
            if os.path.isfile(self.path_system_last):
                with open(self.path_system_last, 'r') as f:
                    self.systemdict.load(f)
        except:
            pass
        
        try:
            if os.path.isfile(self.path_plot_last):
                with open(self.path_plot_last, 'r') as f:
                    self.plotdict.load(f)
        except:
            pass

        # Set validators
        validator_double = DoubleValidator()
        validator_doublenone = DoublenoneValidator()
        validator_doublepositive = DoublepositiveValidator()
        validator_doubledelta = DoubledeltaValidator()
        
        self.ui.lineedit_system_deltaESingle.setValidator(validator_doubledelta)
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
        self.ui.lineedit_system_precision.setValidator(validator_doublepositive)
        self.ui.lineedit_plot_minE.setValidator(validator_doublenone)
        self.ui.lineedit_plot_maxE.setValidator(validator_doublenone)
        
        # Connect signals and slots        
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
        
        self.ui.combobox_system_species1.currentIndexChanged.connect(self.forbidSamebasis)
        self.ui.combobox_system_species2.currentIndexChanged.connect(self.forbidSamebasis)
        
        self.ui.radiobutton_system_pairbasisDefined.toggled.connect(self.togglePairbasis)
        self.ui.radiobutton_plot_overlapDefined.toggled.connect(self.toggleOverlapstate)
        
        self.ui.spinbox_system_deltaNSingle.valueChanged.connect(self.adjustPairlimits)
        self.ui.spinbox_system_deltaLSingle.valueChanged.connect(self.adjustPairlimits)
        self.ui.spinbox_system_deltaJSingle.valueChanged.connect(self.adjustPairlimits)
        self.ui.spinbox_system_deltaMSingle.valueChanged.connect(self.adjustPairlimits)
        
        self.ui.action_system_open.triggered.connect(self.openSystemConf)
        self.ui.action_system_save.triggered.connect(self.saveSystemConf)
        self.ui.action_plot_open.triggered.connect(self.openPlotConf)
        self.ui.action_plot_save.triggered.connect(self.savePlotConf)
        self.ui.action_quit.triggered.connect(self.close)
        self.ui.action_whatsthis.triggered.connect(QtGui.QWhatsThis.enterWhatsThisMode)
        
        self.ui.pushbutton_field1_calc.clicked.connect(self.startCalc)
        self.ui.pushbutton_field2_calc.clicked.connect(self.startCalc)
        self.ui.pushbutton_potential_calc.clicked.connect(self.startCalc)
        
        self.timer.timeout.connect(self.checkForData)
        
        # Emit change-signals in order to let the validation run
        self.ui.spinbox_system_n1.valueChanged.emit(self.ui.spinbox_system_n1.value())
        self.ui.spinbox_system_n2.valueChanged.emit(self.ui.spinbox_system_n2.value())
        self.ui.spinbox_plot_n1.valueChanged.emit(self.ui.spinbox_plot_n1.value())
        self.ui.spinbox_plot_n2.valueChanged.emit(self.ui.spinbox_plot_n2.value())
        
        self.ui.spinbox_system_j1.editingFinished.emit()
        self.ui.spinbox_system_j2.editingFinished.emit()
        self.ui.spinbox_system_m1.editingFinished.emit()
        self.ui.spinbox_system_m2.editingFinished.emit()
        self.ui.spinbox_plot_j1.editingFinished.emit()
        self.ui.spinbox_plot_j2.editingFinished.emit()
        self.ui.spinbox_plot_m1.editingFinished.emit()
        self.ui.spinbox_plot_m2.editingFinished.emit()
        
        self.ui.combobox_system_species1.currentIndexChanged.emit(self.ui.combobox_system_species1.currentIndex())
        
        self.ui.radiobutton_system_pairbasisDefined.toggled.emit(self.ui.radiobutton_system_pairbasisDefined.isChecked())
        self.ui.radiobutton_plot_overlapDefined.toggled.emit(self.ui.radiobutton_plot_overlapDefined.isChecked())
        
        self.ui.spinbox_system_deltaNSingle.valueChanged.emit(self.ui.spinbox_system_deltaNSingle.value())
        self.ui.spinbox_system_deltaLSingle.valueChanged.emit(self.ui.spinbox_system_deltaLSingle.value())
        self.ui.spinbox_system_deltaJSingle.valueChanged.emit(self.ui.spinbox_system_deltaJSingle.value())
        self.ui.spinbox_system_deltaMSingle.valueChanged.emit(self.ui.spinbox_system_deltaMSingle.value())
        
        # Setup plot
        self.minE = None
        self.maxE = None
        
        self.constDistance = self.getConstDistance()
        self.constEField = self.getConstEField()
        self.constBField = self.getConstBField()
        
        for plotarea in [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot, self.ui.graphicsview_potential_plot]:
            plotarea.setDownsampling(ds=True, auto=True, mode='peak')
            plotarea.setClipToView(True)
            plotarea.setLabel('left', 'Energy ('+str(Units.energy)+')')
        
        if self.constEField and not self.constBField:
            for plotarea in [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot]: plotarea.setLabel('bottom', 'Magnetic field ('+str(Units.bfield)+')')
        else:
            for plotarea in [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot]: plotarea.setLabel('bottom', 'Electric field ('+str(Units.efield)+')')
        
        if self.constDistance and not self.constEField:
            self.ui.graphicsview_potential_plot.setLabel('bottom', 'Electric field ('+str(Units.efield)+')')
        elif self.constDistance and not self.constBField:
            self.ui.graphicsview_potential_plot.setLabel('bottom', 'Magnetic field ('+str(Units.bfield)+')')
        else:
            self.ui.graphicsview_potential_plot.setLabel('bottom', 'Interatomic distance ('+str(Units.length)+')')
            
    def resizeEvent(self, event):
        super().resizeEvent(event)
        if sys.platform == "darwin": QtGui.QApplication.processEvents() # hack to circumvent the no-redraw-after-resizing-bug
    
    def getConstEField(self):
        minVec = np.array([self.systemdict['minEx'].magnitude,self.systemdict['minEy'].magnitude,self.systemdict['minEz'].magnitude])
        maxVec = np.array([self.systemdict['maxEx'].magnitude,self.systemdict['maxEy'].magnitude,self.systemdict['maxEz'].magnitude])
        return np.all(minVec==maxVec)
    
    def getConstBField(self):
        minVec = np.array([self.systemdict['minBx'].magnitude,self.systemdict['minBy'].magnitude,self.systemdict['minBz'].magnitude])
        maxVec = np.array([self.systemdict['maxBx'].magnitude,self.systemdict['maxBy'].magnitude,self.systemdict['maxBz'].magnitude])
        return np.all(minVec==maxVec)
    
    def getConstDistance(self):
        minR = self.systemdict['minR'].magnitude # TODO
        maxR = self.systemdict['maxR'].magnitude # TODO
        return minR == maxR
        
    def abortCalculation(self):
        # kill c++ process - this terminates the self.thread, too
        if self.proc is not None: self.proc.kill()
        
        # wait until self.thread has finished
        self.thread.wait()
        
        # clear queues
        self.thread.clear()
    
    def checkForData(self):
        # === print status ===
        elapsedtime = "{}".format(timedelta(seconds=int(time()-self.starttime)))
        if self.thread.message != "":
            self.ui.statusbar.showMessage(self.thread.message+", elapsed time "+elapsedtime)
        else:
            self.ui.statusbar.showMessage("Elapsed time "+elapsedtime)
    
        # === check if memory consumption is to high ===
        if psutil.virtual_memory().percent > 99: # TODO: is the virtual or swap memory the problem on rqo-donkey?
            self.abortCalculation()
            QtGui.QMessageBox.critical(self, "Message", "The program has run out of memory.")
        
        # === process potential map ===
        
        # --- load basis states ---
        
        if self.thread.basisfile_potential != "":
            # load basis
            basis = np.loadtxt(self.thread.basisfile_potential)
                        
            # extract labels
            nlj = basis[:,[1,2,3,5,6,7]]
            
            sorter = np.lexsort(nlj.T[::-1])
            nlj = nlj[sorter]
            diff = np.append([True],np.diff(nlj, axis=0).any(axis=1))
            cumsum = np.cumsum(diff)[np.argsort(sorter)]
            
            self.labelstates_potential = nlj[diff]            
            self.labelmat_potential = sparse.coo_matrix((np.ones_like(cumsum),(np.arange(len(cumsum)),cumsum-1)),shape=(len(cumsum), len(self.labelstates_potential))).tocsc() #nStates, nLabels
            
            self.momentumstrings_potential = [" {}".format(i) for i in np.arange(np.max(self.labelstates_potential[:,[1,4]])+1).astype(np.int)]
            self.momentumstrings_potential[:4] = ['S','P','D','F']
                        
            # remove basis file from hard disk
            os.remove(self.thread.basisfile_potential)
                        
            # indicate that the basis file is already processed
            self.thread.basisfile_potential = ""

        # --- check if there is some new data and if yes, plot it ---
        
        if not self.thread.dataqueue_potential.empty():
        
            # --- storage that allows to draw the hole buffer at once, at least if it is no very large ---
            x = np.array([])
            y = np.array([])
            
            while not self.thread.dataqueue_potential.empty():
                # --- load eigenvalues (energies, y value) and eigenvectors (basis) ---
                filestep, filename = self.thread.dataqueue_potential.get()
                eigensystem = Eigensystem(filename)
                energies = eigensystem.energies
                basis = eigensystem.basis
            
                # --- determine which basis elements are within the energy range ---
                boolarr = np.ones(len(energies),dtype=np.bool)
                if self.minE is not None: boolarr &= energies >= self.minE
                if self.maxE is not None: boolarr &= energies <= self.maxE
            
                # cut the energies
                energies = energies[boolarr]
            
                # convert the energies
                energies *= self.converter_y
            
                # cut the basis
                idxarr, = np.nonzero(boolarr)
                transformator = sparse.coo_matrix((np.ones_like(idxarr),(idxarr,np.arange(len(idxarr)))),shape=(basis.shape[1], len(idxarr))).tocsc()
                basis *= transformator
                
                # probability to be in a certain state
                probs = np.abs(basis).power(2) # nState, nBasis
            
                # --- calculate the position (x value) ---
                if self.constDistance and not self.constEField:
                    vec = np.array([float(eigensystem.params["Ex"]),float(eigensystem.params["Ey"]),float(eigensystem.params["Ez"])])
                    position = np.sign(np.vdot(vec,[1,1,1]))*np.linalg.norm(vec)*self.converter_x_potential
                elif self.constDistance and not self.constBField:
                    vec = np.array([float(eigensystem.params["Bx"]),float(eigensystem.params["By"]),float(eigensystem.params["Bz"])])
                    position = np.sign(np.vdot(vec,[1,1,1]))*np.linalg.norm(vec)*self.converter_x_potential
                else:
                    position = float(eigensystem.params["R"])*self.converter_x_potential
            
                # --- store data to plot several points at once ---  
                x = np.append(x,position*np.ones_like(energies))
                y = np.append(y,energies)
            
                # --- draw labels at the beginning of the plotting --- 
                if self.ui.groupbox_plot_labels.isChecked() and filestep == 0:
            
                    # probability to find a label inside a basis element
                    labelprob = probs.T*self.labelmat_potential #nBasis, nLabels
                
                    # total probability to find a label
                    cumprob = np.array(labelprob.sum(axis=0).flat)
                    boolarr = cumprob > 0.1
                
                    # normalize in such a way that the total probability is always one
                    idxarr, = np.nonzero(boolarr)
                    normalizer = sparse.coo_matrix((1/cumprob[idxarr],(idxarr,np.arange(len(idxarr)))),shape=(labelprob.shape[1], len(idxarr))).tocsc()
                
                    # store the energetic expectation value of the labels
                    labelenergies = ((labelprob*normalizer).T*energies)
                
                    # store the position of the labels
                    labelposition = position

                    # get size and alpha value of labels
                    size = '{}'.format(max(int(round(self.ui.spinbox_plot_szLabel.value()*11)),1))
                    alpha = int(round(self.ui.spinbox_plot_transpLabel.value()*255))
                
                    # draw the labels
                    for labelstate, labelenergy in zip(self.labelstates_potential[boolarr],labelenergies):
                        sn1, sl1, sj1, sn2, sl2, sj2 = labelstate

                        text = pg.TextItem(html='<div style="text-align: center; font-size: '+size+'pt;"><span style="color: rgba(0,0,0,255);">'\
                            +'{}{}<sub style="font-size: '.format(int(sn1),self.momentumstrings_potential[int(sl1)]) \
                            +size+'pt;">{}/2</sub>'.format(int(2*sj1))\
                            +' {}{}<sub style="font-size: '.format(int(sn2),self.momentumstrings_potential[int(sl2)]) \
                            +size+'pt;">{}/2</sub>'.format(int(2*sj2))\
                            +'</span></div>',anchor=(1, 0.5),fill=(250,235,215,alpha),border=(255,228,181,255)) # TODO anchor dependent on start < end !!!!!!!!!!!!!!!!!!!
                        text.setPos(labelposition, labelenergy)
                        text.setZValue(15)
                        self.ui.graphicsview_potential_plot.addItem(text)

                    curve = PointsItem(labelposition*np.ones_like(labelenergies), labelenergies, 0, 0, (0,0,0))
                    curve.setZValue(5)
                    self.ui.graphicsview_potential_plot.addItem(curve)
                
                    # break since drawing labels already take some time
                    break
            
                # --- break if enough points are collected --- 
                if len(x) > 5000: break                
        
            # --- plot the stored data ---
            if self.ui.groupbox_plot_points.isChecked() and len(x) > 0:
        
                # get size and alpha value of points
                size = self.ui.spinbox_plot_szPoint.value()
                alpha = self.ui.spinbox_plot_transpPoint.value()*255
                
                # plot the basis elements
                curve = PointsItem(x, y, size, alpha, (0,0,0))
                curve.setZValue(10)
                self.ui.graphicsview_potential_plot.addItem(curve)
        
            # --- update the graphics view ---
            self.ui.graphicsview_potential_plot.repaint()
        
        # === process field maps ===
             
        for idx in range(2):
        
            dataqueue = [self.thread.dataqueue_field1, self.thread.dataqueue_field2][idx]
            basisfile = [self.thread.basisfile_field1, self.thread.basisfile_field2][idx]
            graphicsview_plot = [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot]
            
            # --- load basis states ---
            
            if basisfile != "":
                # load basis
                basis = np.loadtxt(basisfile)
            
                # calculate a matrix that can be used to determine the momenta inside a basis element 
                momentum = basis[:,2]
                self.momentummat[idx] = sparse.csc_matrix((momentum[:,None] == np.arange(np.max(momentum)+1)[None,:]).astype(int))

                # extract labels
                nlj = basis[:,[1,2,3]]
                
                sorter = np.lexsort(nlj.T[::-1])
                nlj = nlj[sorter]
                diff = np.append([True],np.diff(nlj, axis=0).any(axis=1))
                cumsum = np.cumsum(diff)[np.argsort(sorter)]
            
                self.labelstates[idx] = nlj[diff]            
                self.labelmat[idx] = sparse.coo_matrix((np.ones_like(cumsum),(np.arange(len(cumsum)),cumsum-1)),shape=(len(cumsum), len(self.labelstates[idx]))).tocsc() #nStates, nLabels
                #self.labelstates[idx] = unique_rows(nlj)
                #self.labelmat[idx] = sparse.csc_matrix((nlj[:,None,:] == self.labelstates[idx][None,:,:]).all(axis=-1).astype(int)) #nStates, nLabels
            
                self.momentumstrings[idx] = [" {}".format(i) for i in np.arange(np.max(self.labelstates[idx][:,1])+1).astype(np.int)]
                self.momentumstrings[idx][:4] = ['S','P','D','F']  
            
                # remove basis file from hard disk
                os.remove(basisfile)
            
                # indicate that the basis file is already processed
                if idx == 0: self.thread.basisfile_field1 = ""
                elif idx == 1: self.thread.basisfile_field2 = ""
                
            
            # --- check if there is some new data and if yes, plot it ---
            
            if not dataqueue.empty():
            
                # --- storage that allows to draw the hole buffer at once, at least if it is no very large ---
                x = np.array([])
                y = np.array([])
                l = np.array([])
            
                while not dataqueue.empty():
                    # --- load eigenvalues (energies, y value) and eigenvectors (basis) ---
                    filestep, filename = dataqueue.get()
                    eigensystem = Eigensystem(filename)
                    energies = eigensystem.energies
                    basis = eigensystem.basis
                
                    # --- determine which basis elements are within the energy range ---
                    boolarr = np.ones(len(energies),dtype=np.bool)
                    if self.minE is not None: boolarr &= energies >= self.minE
                    if self.maxE is not None: boolarr &= energies <= self.maxE
                
                    # cut the energies
                    energies = energies[boolarr]
                
                    # convert the energies
                    energies *= self.converter_y
                
                    # cut the basis
                    idxarr, = np.nonzero(boolarr)
                    transformator = sparse.coo_matrix((np.ones_like(idxarr),(idxarr,np.arange(len(idxarr)))),shape=(basis.shape[1], len(idxarr))).tocsc()
                    basis *= transformator
                    
                    # probability to be in a certain state
                    probs = np.abs(basis).power(2) # nState, nBasis
                    
                    # --- calculate the momentum that is associated with a basis element ---
                    # calculate which momenta appear in a basis element
                    momentum_probabilty = probs.T*self.momentummat[idx] # nBasis, nMomentum
                
                    # keep information about the momentum that appears with a probability > 0.5 only
                    momentum_probabilty.data[:] *= momentum_probabilty.data>0.5
                    momentum_probabilty.eliminate_zeros()
                    momentum_probabilty = sparse.coo_matrix(momentum_probabilty)
                
                    # store the value of this momentum
                    momentum = -np.ones(momentum_probabilty.shape[0],dtype=np.int) # -1 means no determinable momentum
                    momentum[momentum_probabilty.row] = momentum_probabilty.col
                
                    # --- calculate the position (x value) ---                
                    if self.constEField and not self.constBField:
                        vec = np.array([float(eigensystem.params["Bx"]),float(eigensystem.params["By"]),float(eigensystem.params["Bz"])])
                        position = np.sign(np.vdot(vec,[1,1,1]))*np.linalg.norm(vec)*self.converter_x_field
                    else:
                        vec = np.array([float(eigensystem.params["Ex"]),float(eigensystem.params["Ey"]),float(eigensystem.params["Ez"])])
                        position = np.sign(np.vdot(vec,[1,1,1]))*np.linalg.norm(vec)*self.converter_x_field
                
                    # --- store data to plot several points at once ---  
                    x = np.append(x,position*np.ones_like(energies))
                    y = np.append(y,energies)
                    l = np.append(l,momentum)
                
                    # --- draw labels at the beginning of the plotting --- 
                    if self.ui.groupbox_plot_labels.isChecked() and filestep == 0:
                
                        # probability to find a label inside a basis element
                        labelprob = probs.T*self.labelmat[idx] #nBasis, nLabels
                    
                        # total probability to find a label
                        cumprob = np.array(labelprob.sum(axis=0).flat)
                        boolarr = cumprob > 0.1
                    
                        # normalize in such a way that the total probability is always one
                        idxarr, = np.nonzero(boolarr)
                        normalizer = sparse.coo_matrix((1/cumprob[idxarr],(idxarr,np.arange(len(idxarr)))),shape=(labelprob.shape[1], len(idxarr))).tocsc()
                    
                        # store the energetic expectation value of the labels
                        labelenergies = ((labelprob*normalizer).T*energies)
                    
                        # store the position of the labels
                        labelposition = position

                        # get size and alpha value of labels
                        size = '{}'.format(max(int(round(self.ui.spinbox_plot_szLabel.value()*11)),1))
                        alpha = int(round(self.ui.spinbox_plot_transpLabel.value()*255))
                    
                        # draw the labels
                        for labelstate, labelenergy in zip(self.labelstates[idx][boolarr],labelenergies):
                            sn, sl, sj = labelstate
                            
                            for j in range(2):
                                if not self.thread.samebasis and j != idx: continue
                                text = pg.TextItem(html='<div style="text-align: center; font-size: '+size+'pt;">'\
                                    +'<span style="color: rgba(0,0,0,255);">{}{}<sub style="font-size: '.format(int(sn),self.momentumstrings[idx][int(sl)]) \
                                    +size+'pt;">{}/2</sub></span></div>'.format(int(2*sj)),anchor=(0, 0.5),fill=(250,235,215,alpha),border=(255,228,181,255)) # TODO anchor dependent on start < end !!!!!!!!!!!!!!!!!!!
                                text.setPos(labelposition, labelenergy)
                                text.setZValue(15)
                                graphicsview_plot[j].addItem(text)
                            
                        for j in range(2):
                            if not self.thread.samebasis and j != idx: continue
                            curve = PointsItem(labelposition*np.ones_like(labelenergies), labelenergies, 0, 0, (0,0,0))
                            curve.setZValue(5)
                            graphicsview_plot[j].addItem(curve)
                    
                        # break since drawing labels already take some time
                        break
                
                    # --- break if enough points are collected --- 
                    if len(x) > 5000: break                
            
                # --- plot the stored data ---
                if self.ui.groupbox_plot_points.isChecked() and len(x) > 0:
            
                    # get size and alpha value of points
                    size = self.ui.spinbox_plot_szPoint.value()
                    alpha = self.ui.spinbox_plot_transpPoint.value()*255
                
                    # cut the momenta to reasonable values
                    l[l>len(self.momentumcolors)-2] = len(self.momentumcolors)-2
                    l[l<0] = len(self.momentumcolors)-1
                
                    # loop over momenta
                    for i in range(len(self.momentumcolors)):
                        # determine which basis elements have the current momentum
                        boolarr = l == i
                        if (np.sum(boolarr) == 0): continue
                    
                        # determine the associated color
                        color = self.momentumcolors[i]
                    
                        # plot the basis elements
                        for j in range(2):
                            if not self.thread.samebasis and j != idx: continue
                            curve = PointsItem(x[boolarr], y[boolarr], size, alpha, color)
                            curve.setZValue(10)
                            graphicsview_plot[j].addItem(curve)
            
                # --- update the graphics view ---
                for j in range(2):
                    if not self.thread.samebasis and j != idx: continue
                    graphicsview_plot[j].repaint()
        

            
            
            
            
            
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
            
        
        # === terminate thread ===
        
        # check if thread has finished
        if self.thread.isFinished() and self.thread.dataqueue_field1.empty() and self.thread.dataqueue_field2.empty() and self.thread.dataqueue_potential.empty():
            
            """c6 = 0.008750183201815791 # 40
            c6 = 0.009304149469320561 # 39
            c6 = 0.009104493144885473 # 38
            c6 = -0.044880307113037206 # 46
            c6 = -465453.7925703782

            #c6 = -10.306017051622977
            #c6 = -450153.7925703782
            x = np.linspace(self.systemdict["minR"].magnitude,self.systemdict["maxR"].magnitude,500)
            y = -c6/x**6
            
            if self.plotdict["minE"] is not None:
                minE = self.plotdict["minE"].magnitude
                x = x[y >= minE]
                y = y[y >= minE]
            if self.plotdict["maxE"] is not None:
                maxE = self.plotdict["maxE"].magnitude
                x = x[y <= maxE]
                y = y[y <= maxE]
            
            self.ui.graphicsview_potential_plot.plot(x,y,pen={'color': (255,0,0), 'width': 1})"""
            
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
            if self.ui.pushbutton_field1_calc == self.senderbutton:
                self.senderbutton.setText("Calculate field map of atom 1")
            if self.ui.pushbutton_field2_calc == self.senderbutton:
                self.senderbutton.setText("Calculate field map of atom 2")
            if self.ui.pushbutton_potential_calc == self.senderbutton:
                self.senderbutton.setText("Calculate pair potential")

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
    
    @QtCore.pyqtSlot(float)
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
    
    @QtCore.pyqtSlot(bool)
    def togglePairbasis(self):
        checked = self.ui.radiobutton_system_pairbasisDefined.isChecked()
        self.ui.widget_system_pair.setEnabled(checked)
    
    @QtCore.pyqtSlot(bool)
    def toggleOverlapstate(self):
        checked = self.ui.radiobutton_plot_overlapDefined.isChecked()
        self.ui.widget_plot_qn.setEnabled(checked)
    
    @QtCore.pyqtSlot(str)
    def forbidSamebasis(self):
        if self.ui.combobox_system_species1.currentIndex() != self.ui.combobox_system_species2.currentIndex():
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
                
                if n_err or l_err or j_err or m_err:
                    self.invalidQuantumnumbers[i] = True
                else:
                    self.invalidQuantumnumbers[i] = False
        
        if np.any(self.invalidQuantumnumbers):
            self.ui.statusbar.showMessage('Invalide quantum numbers specified.')
        else:
            self.ui.statusbar.showMessage('')
    
    @QtCore.pyqtSlot(str)   
    def validateHalfinteger(self):
        value = self.sender().value()
        self.sender().setValue(np.floor(value)+0.5)
            
    @QtCore.pyqtSlot()
    def startCalc(self):         
        if self.proc is None:
            if np.any(self.invalidQuantumnumbers):
                QtGui.QMessageBox.critical(self, "Message", "Invalide quantum numbers specified.")
                
            else:
                # Save last settings
                with open(self.path_system_last, 'w') as f:
                    self.systemdict.saveInOriginalunits(f)
        
                with open(self.path_plot_last, 'w') as f:
                    self.plotdict.saveInOriginalunits(f)
            
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
                self.constDistance = self.getConstDistance()
                self.constEField = self.getConstEField()
                self.constBField = self.getConstBField()
            
                if self.senderbutton in [self.ui.pushbutton_field1_calc, self.ui.pushbutton_potential_calc] or self.samebasis:
                    self.ui.graphicsview_field1_plot.clear()
                    self.ui.graphicsview_field1_plot.setLabel('left', 'Energy ('+str(Units.energy)+')')
                    if self.constEField and not self.constBField:
                        self.ui.graphicsview_field1_plot.setLabel('bottom', 'Magnetic field ('+str(Units.bfield)+')')
                        self.converter_x_field = Converter.fromAU(1,Units.bfield).magnitude
                    else:
                        self.ui.graphicsview_field1_plot.setLabel('bottom', 'Electric field ('+str(Units.efield)+')')
                        self.converter_x_field = Converter.fromAU(1,Units.efield).magnitude
                if self.senderbutton in [self.ui.pushbutton_field2_calc, self.ui.pushbutton_potential_calc] or self.samebasis:
                    self.ui.graphicsview_field2_plot.clear()
                    self.ui.graphicsview_field2_plot.setLabel('left', 'Energy ('+str(Units.energy)+')')
                    if self.constEField and not self.constBField:
                        self.ui.graphicsview_field2_plot.setLabel('bottom', 'Magnetic field ('+str(Units.bfield)+')')
                        self.converter_x_field = Converter.fromAU(1,Units.bfield).magnitude
                    else:
                        self.ui.graphicsview_field2_plot.setLabel('bottom', 'Electric field ('+str(Units.efield)+')')
                        self.converter_x_field = Converter.fromAU(1,Units.efield).magnitude
                if self.senderbutton == self.ui.pushbutton_potential_calc:
                    self.ui.graphicsview_potential_plot.clear()
                    self.ui.graphicsview_potential_plot.setLabel('left', 'Energy ('+str(Units.energy)+')')
                    if self.constDistance and not self.constEField:
                        self.ui.graphicsview_potential_plot.setLabel('bottom', 'Electric field ('+str(Units.efield)+')')
                        self.converter_x_potential = Converter.fromAU(1,Units.efield).magnitude
                    elif self.constDistance and not self.constBField:
                        self.ui.graphicsview_potential_plot.setLabel('bottom', 'Magnetic field ('+str(Units.bfield)+')')
                        self.converter_x_potential = Converter.fromAU(1,Units.bfield).magnitude
                    else:
                        self.ui.graphicsview_potential_plot.setLabel('bottom', 'Interatomic distance ('+str(Units.length)+')')
                        self.converter_x_potential = Converter.fromAU(1,Units.length).magnitude
                        
                self.converter_y = Converter.fromAU(1,Units.energy).magnitude 
                    
                # Set limits
                self.minE = self.plotdict["minE"]
                if self.minE is not None:
                    self.minE = Converter.toAU(self.minE).magnitude
                
                self.maxE = self.plotdict["maxE"]
                if self.maxE is not None:
                    self.maxE = Converter.toAU(self.maxE).magnitude
            
                # Save configuration to json file
                with open(self.path_config, 'w') as f:
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
                    
                # Start c++ process
                if self.systemdict["minEx"].magnitude != 0 or self.systemdict["minEy"].magnitude != 0 or self.systemdict["maxEx"].magnitude != 0 or self.systemdict["maxEy"].magnitude != 0 or \
                        self.systemdict["minBx"].magnitude != 0 or self.systemdict["minBy"].magnitude != 0 or self.systemdict["maxBx"].magnitude != 0 or self.systemdict["maxBy"].magnitude != 0:
                    path_cpp = self.path_cpp_complex
                else:
                    path_cpp = self.path_cpp_real
                    
                self.proc = subprocess.Popen(["mpiexec","-n","%d"%self.numprocessors,path_cpp,"-c",self.path_config, "-o", self.path_out],
                    stdout=subprocess.PIPE, cwd=self.path_workingdir)
                
                self.starttime = time()
        
                # Start thread that collects the output
                self.thread.execute(self.proc.stdout)
            
                # Start timer used for processing the results
                self.timer.start(0)
            
        else:
            self.abortCalculation()            
    
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
        self.abortCalculation()
        
        # Save last settings
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
