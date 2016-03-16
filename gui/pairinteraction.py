#!/usr/bin/env python3



import sys
from pint import UnitRegistry
from pint.unit import UndefinedUnitError
from PyQt4 import QtCore, QtGui
from plotter import Ui_plotwindow # pyuic4 plotter.ui > plotter.py or py3uic4 plotter.ui > plotter.py
import pyqtgraph as pg
import numpy as np
import collections
from abc import ABCMeta, abstractmethod
from time import sleep, time
from datetime import timedelta
import locale
import json
import os
import multiprocessing
import subprocess
from scipy import sparse
import signal
from queue import Queue
import ctypes
import sip
from scipy import constants
import psutil
from scipy.ndimage.filters import gaussian_filter
from palettable import cubehelix

locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
QtCore.QLocale.setDefault(QtCore.QLocale(locale.getlocale()[0]))

signal.signal(signal.SIGINT, signal.SIG_DFL)

ureg = UnitRegistry()
Q = ureg.Quantity
C = lambda s : Q(constants.value(s),constants.unit(s))

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

# http://stackoverflow.com/questions/4695337/expanding-adding-a-row-or-column-a-scipy-sparse-matrix
def csr_vappend(a,b):
    """ Takes in 2 csr_matrices and appends the second one to the bottom of the first one. 
    Much faster than scipy.sparse.vstack but assumes the type to be csr and overwrites
    the first matrix instead of copying it. The data, indices, and indptr still get copied."""

    a.data = np.hstack((a.data,b.data))
    a.indices = np.hstack((a.indices,b.indices))
    a.indptr = np.hstack((a.indptr,(b.indptr + a.nnz)[1:]))
    a._shape = (a.shape[0]+b.shape[0],b.shape[1])
    #a = sparse.csr_matrix( (a.data,a.indices,a.indptr), shape=a.shape )


# ----- Returns a normalized image -----

def normscale(data, cmin=None, cmax=None):
    if cmin is None:
        cmin = np.nanmin(data)
    if cmax is None:
        cmax = np.nanmax(data)
    return (data-cmin)/(cmax-cmin or 1)

# ----- Returns a byte-scaled image -----

def bytescale(data, cmin=None, cmax=None, high=255, low=0):
    return np.array(normscale(data,cmin,cmax)*(high - low + 0.9999) + low).astype(int) # TODO


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
    
    def load(self, params):
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
    
    def paramsInOriginalunits(self, exclude = []):
        params = dict()
        for k, v in self.items():
            if k in exclude: continue
            if isinstance(v, Q): params[k] = str(v)
            else: params[k] = v
        return params

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
        store["minE_field1"] = {'widget': ui.lineedit_field1_minE, 'unit': Units.energy}
        store["maxE_field1"] = {'widget': ui.lineedit_field1_maxE, 'unit': Units.energy}
        store["minE_field2"] = {'widget': ui.lineedit_field2_minE, 'unit': Units.energy}
        store["maxE_field2"] = {'widget': ui.lineedit_field2_maxE, 'unit': Units.energy}
        store["minE_potential"] = {'widget': ui.lineedit_potential_minE, 'unit': Units.energy}
        store["maxE_potential"] = {'widget': ui.lineedit_potential_maxE, 'unit': Units.energy}
        store["lines"] = {'widget': ui.groupbox_plot_lines}
        store["points"] = {'widget': ui.groupbox_plot_points}
        store["labels"] = {'widget': ui.groupbox_plot_labels}
        store["overlap"] = {'widget': ui.groupbox_plot_overlap}
        store["szLine"] = {'widget': ui.spinbox_plot_szLine, 'unit': Units.dimensionless}
        store["szPoint"] = {'widget': ui.spinbox_plot_szPoint, 'unit': Units.dimensionless}
        store["szLabel"] = {'widget': ui.spinbox_plot_szLabel, 'unit': Units.dimensionless}
        store["szOverlap"] = {'widget': ui.spinbox_plot_szOverlap, 'unit': Units.dimensionless}
        store["transpLine"] = {'widget': ui.spinbox_plot_transpLine, 'unit': Units.dimensionless}
        store["transpPoint"] = {'widget': ui.spinbox_plot_transpPoint, 'unit': Units.dimensionless}
        store["transpLabel"] = {'widget': ui.spinbox_plot_transpLabel, 'unit': Units.dimensionless}
        store["transpOverlap"] = {'widget': ui.spinbox_plot_transpOverlap, 'unit': Units.dimensionless}
        store["overlapUnperturbed"] = {'widget': ui.radiobutton_plot_overlapUnperturbed}
        store["lin"] = {'widget': ui.radiobutton_plot_lin}
        store["log"] = {'widget': ui.radiobutton_plot_log}
        store["resolution"] = {'widget': ui.spinbox_plot_resolution}
        store["n1"] = {'widget': ui.spinbox_plot_n1, 'unit': Units.dimensionless}
        store["n2"] = {'widget': ui.spinbox_plot_n2, 'unit': Units.dimensionless}
        store["l1"] = {'widget': ui.spinbox_plot_l1, 'unit': Units.dimensionless}
        store["l2"] = {'widget': ui.spinbox_plot_l2, 'unit': Units.dimensionless}
        store["j1"] = {'widget': ui.spinbox_plot_j1, 'unit': Units.dimensionless}
        store["j2"] = {'widget': ui.spinbox_plot_j2, 'unit': Units.dimensionless}
        store["m1"] = {'widget': ui.spinbox_plot_m1, 'unit': Units.dimensionless}
        store["m2"] = {'widget': ui.spinbox_plot_m2, 'unit': Units.dimensionless}
        store["antialiasing"] = {'widget': ui.checkbox_plot_antialiasing}

class Worker(QtCore.QThread):

    def __init__(self, parent = None):
        super().__init__(parent)
        self.exiting = False
        self.samebasis = False
        self.message = ""
        self.basisfile_field1 = ""
        self.basisfile_field2 = ""
        self.basisfile_potential = ""
        self.numBlocks_field1 = 0
        self.numBlocks_field2 = 0
        self.numBlocks_potential = 0
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
                numBlocks = int(line[12:19].decode('utf-8'))
                current = 0
                if type == 0 or type == 3:
                    self.numBlocks_field1 = numBlocks
                elif type == 1:
                    self.numBlocks_field2 = numBlocks
                elif type == 2:
                    self.numBlocks_potential = numBlocks
                
            elif line[:5] == b">>DIM":
                dim = int(line[5:12])
                status_progress = "diagonalize {} x {} matrix, {} of {} matrices processed".format(dim, dim, current,total)
                
            elif line[:5] == b">>OUT":
                current += 1
                status_progress = "diagonalize {} x {} matrix, {} of {} matrices processed".format(dim, dim, current,total)
                
                filenumber = int(line[5:12].decode('utf-8'))
                filestep = int(line[12:19].decode('utf-8'))
                blocknumber = int(line[19:26].decode('utf-8'))
                filename = line[27:-1].decode('utf-8')
                
                if type == 0 or type == 3:
                    self.dataqueue_field1.put([filestep,blocknumber,filename])
                elif type == 1:
                    self.dataqueue_field2.put([filestep,blocknumber,filename])
                elif type == 2:
                    self.dataqueue_potential.put([filestep,blocknumber,filename])
                    
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
    def __init__(self, x, y, size=1, alpha=80, color=(0,0,0)):
        """x and y are 2D arrays of shape (Nplots, Nsamples)"""
        connections = np.ones(x.shape, dtype=bool)
        connections[:,-1] = 0 # don't draw the segment between each trace
        
        self.path = pg.arrayToQPath(x.flatten(), y.flatten(), connections.flatten())
        pg.QtGui.QGraphicsPathItem.__init__(self, self.path)
        pen = pg.mkPen(color+(alpha,),width=size,cosmetic=True)
        self.setPen(pen)
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
        
        if status[0] == QtGui.QValidator.Acceptable and locale.atof(s) < 0: # TODO atof soll auch mit englischem Zahlenformat umgehen koennen
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
        
        del pg.graphicsItems.GradientEditorItem.Gradients['greyclip']
        del pg.graphicsItems.GradientEditorItem.Gradients['grey']
        del pg.graphicsItems.GradientEditorItem.Gradients['cyclic']
        del pg.graphicsItems.GradientEditorItem.Gradients['spectrum']
        del pg.graphicsItems.GradientEditorItem.Gradients['bipolar']
        
        
        color = cubehelix.Cubehelix.make(sat=1.8,n=7,rotation=1.21,start=1.2,reverse=True).colors[::-1]
        color = np.append(color, [[255]]*len(color), axis=1).astype(np.ubyte)
        pos = np.linspace(0,1,len(color))
        pg.graphicsItems.GradientEditorItem.Gradients['cubehelix1'] = {'mode':'rgb', 'ticks':[(p, tuple(c)) for c, p in zip(color,pos)]}
        
        
        color = cubehelix.Cubehelix.make(sat=1.5,n=7,rotation=-1.0,start=0.9,reverse=True).colors[::-1]
        color = np.append(color, [[255]]*len(color), axis=1).astype(np.ubyte)
        pos = np.linspace(0,1,len(color))
        pg.graphicsItems.GradientEditorItem.Gradients['cubehelix2'] = {'mode':'rgb', 'ticks':[(p, tuple(c)) for c, p in zip(color,pos)]}
        
        for k, v in pg.graphicsItems.GradientEditorItem.Gradients.items():
            pg.graphicsItems.GradientEditorItem.Gradients[k]['ticks'] = [(1-p, c) for p, c in v['ticks']]
        
        
        
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
        self.path_tabs_last = os.path.join(self.path_out,"lasttabs.json")
        self.path_config = os.path.join(self.path_out,"conf.json")
        
        self.proc = None
        
        self.thread = Worker()
        self.timer = QtCore.QTimer()
        
        self.momentumcolors = [(55,126,184),(77,175,74),(228,26,28),(152,78,163),(0,0,0),(255//5,255//5,255//5)] # s, p, d, f, other, undetermined
        self.symmetrycolors = [(0,0,0),(140,81,10),(1,102,94)] # all, sym, asym
        
        self.momentummat = [None]*2
        self.labelmat = [None]*2
        self.labelstates = [None]*2
        self.momentumstrings = [None]*2
        
        self.buffer_basis = [{}]*2
        self.buffer_energies = [{}]*2
        self.buffer_positions = [{}]*2
        self.buffer_boolarr = [{}]*2
        self.buffer_basis_potential = {}
        self.buffer_energies_potential = {}
        self.buffer_positions_potential = {}
        
        self.buffer_energiesMap = [{}]*2
        self.buffer_positionsMap = [{}]*2
        self.buffer_energiesMap_potential = {}
        self.buffer_positionsMap_potential = {}
        self.buffer_overlapMap_potential = {}
        
        self.labelprob_potential = None
        self.labelprob_num_potential = 0
        
        self.lines_buffer_minIdx = {}
        self.colormap_buffer_minIdx = 0
        self.lines_buffer_minIdx_field = 0
        

        
        
        #clrmp = pg.ColorMap(pos,color)
        #self.lut = clrmp.getLookupTable()
                
        self.ui.gradientwidget_plot_gradient.setOrientation("top")
        self.ui.gradientwidget_plot_gradient.loadPreset('cubehelix1')

        # TODOs
        self.ui.lineedit_system_theta.setEnabled(False)
        self.ui.lineedit_system_precision.setEnabled(False)
        self.ui.pushbutton_field1_save.setEnabled(False)
        self.ui.pushbutton_field2_save.setEnabled(False)
        self.ui.pushbutton_potential_save.setEnabled(False)
        self.ui.checkbox_system_dq.setEnabled(False)
        self.ui.checkbox_system_qq.setEnabled(False)
        
        # Create directories
        if not os.path.exists(self.path_out):
            os.makedirs(self.path_out)	
            if os.name == 'nt':
                ret = ctypes.windll.kernel32.SetFileAttributesW(self.path_out,FILE_ATTRIBUTE_HIDDEN)
                if not ret: raise ctypes.WinError()
                
        

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
        self.ui.lineedit_field1_minE.setValidator(validator_doublenone)
        self.ui.lineedit_field1_maxE.setValidator(validator_doublenone)
        self.ui.lineedit_field2_minE.setValidator(validator_doublenone)
        self.ui.lineedit_field2_maxE.setValidator(validator_doublenone)
        self.ui.lineedit_potential_minE.setValidator(validator_doublenone)
        self.ui.lineedit_potential_maxE.setValidator(validator_doublenone)
        
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
        self.ui.radiobutton_plot_log.toggled.connect(self.toggleYScale)
        
        self.ui.checkbox_plot_antialiasing.toggled.connect(self.toggleAntialiasing)
        
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
        
        # Load last settings
        try: # TODO
            if os.path.isfile(self.path_system_last):
                with open(self.path_system_last, 'r') as f:
                    params = json.load(f)
                    self.systemdict.load(params)
        except:
            pass
        
        try:
            if os.path.isfile(self.path_plot_last):
                with open(self.path_plot_last, 'r') as f:
                    params = json.load(f)
                    self.ui.gradientwidget_plot_gradient.restoreState(params["gradientwidget"])
                    del params["gradientwidget"]
                    self.plotdict.load(params)
        except:
            pass
        
        try:
            if os.path.isfile(self.path_tabs_last):
                with open(self.path_tabs_last, 'r') as f:
                    params = json.load(f)
                    self.ui.tabwidget_config.setCurrentIndex(params["config"])
                    self.ui.tabwidget_plotter.setCurrentIndex(params["plotter"])
                    self.ui.toolbox_system.setCurrentIndex(params["system"])
        except:
            pass
        
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
        self.ui.radiobutton_plot_log.toggled.emit(self.ui.radiobutton_plot_log.isChecked())
        
        self.ui.checkbox_plot_antialiasing.toggled.emit(self.ui.checkbox_plot_antialiasing.isChecked())
        
        self.ui.spinbox_system_deltaNSingle.valueChanged.emit(self.ui.spinbox_system_deltaNSingle.value())
        self.ui.spinbox_system_deltaLSingle.valueChanged.emit(self.ui.spinbox_system_deltaLSingle.value())
        self.ui.spinbox_system_deltaJSingle.valueChanged.emit(self.ui.spinbox_system_deltaJSingle.value())
        self.ui.spinbox_system_deltaMSingle.valueChanged.emit(self.ui.spinbox_system_deltaMSingle.value())
        
        
        
        
        # Setup plot
        self.minE_field1 = None
        self.minE_field2 = None
        self.minE_potential = None
        self.maxE_field1 = None
        self.maxE_field2 = None
        self.maxE_potential = None
        
        self.constDistance = self.getConstDistance()
        self.constEField = self.getConstEField()
        self.constBField = self.getConstBField()
        #self.sameSpecies = self.getSameSpecies()
        
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
    
    def get1DPosition(self, x, y, z):
        vec = np.array([x,y,z])
        return np.sign(np.vdot(vec,[1,1,1]))*np.linalg.norm(vec)
        
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
    
    #def getSameSpecies(self):
    #    return self.systemdict['species1'] == self.systemdict['species2']
        
    def abortCalculation(self):
        # kill c++ process - this terminates the self.thread, too
        if self.proc is not None: self.proc.kill()
        
        # wait until self.thread has finished
        self.thread.wait()
        
        # clear queues
        self.thread.clear()
    
    def checkForData(self):
        dataamount = 0
    
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
            
            if self.ui.groupbox_plot_overlap.isChecked():
                self.stateidx_potential = np.where(np.all(basis[:,[1,2,3,4,5,6,7,8]] == self.overlapstate[None,:],axis=-1))[0]
                if len(self.stateidx_potential) == 0: self.stateidx_potential =  -1
                else: self.stateidx_potential = self.stateidx_potential[0]
                
                self.displayunits2pixelunits_x = np.nan
                self.displayunits2pixelunits_y = np.nan
                        
            # extract labels
            nlj = basis[:,[1,2,3,5,6,7]]
            
            # sort pair state names
            if self.thread.samebasis:
                firstsmaller = np.argmax(np.append(nlj[:,0:3] < nlj[:,3:6],np.ones((len(nlj),1),dtype=bool),axis=-1),axis=-1) # TODO in Funktion auslagern
                firstsmaller_reverted = np.argmax(np.append(nlj[:,3:6] < nlj[:,0:3],np.ones((len(nlj),1),dtype=bool),axis=-1),axis=-1) # TODO in Funktion auslagern
                namesToSwap = firstsmaller > firstsmaller_reverted
                nlj[namesToSwap] = nlj[namesToSwap][:,[3,4,5,0,1,2]]
            
            sorter = np.lexsort(nlj.T[::-1])
            nlj = nlj[sorter]
            
            diff = np.append([True],np.diff(nlj, axis=0).any(axis=1))
            
            """if self.sameSpecies:
                order = np.sum(nlj * (np.array([1,2,3,4,5,6]) * 10000)[None,:], axis = -1)
                order_reversed = np.sum(nlj * (np.array([4,5,6,1,2,3]) * 10000)[None,:], axis = -1)
                print(np.sum(diff), np.sum(order <= order_reversed), len(order))
                diff &= order <= order_reversed
                print(np.sum(diff), np.sum(order <= order_reversed), len(order))"""
            
            cumsum = np.cumsum(diff)[np.argsort(sorter)]
            
            self.labelstates_potential = nlj[diff]            
            self.labelmat_potential = sparse.coo_matrix((np.ones_like(cumsum),(np.arange(len(cumsum)),cumsum-1)),shape=(len(cumsum), len(self.labelstates_potential))).tocsr() #nStates, nLabels
                
            self.momentumstrings_potential = [" {}".format(i) for i in np.arange(np.max(self.labelstates_potential[:,[1,4]])+1).astype(np.int)]
            self.momentumstrings_potential[:4] = ['S','P','D','F']
                        
            # remove basis file from hard disk
            os.remove(self.thread.basisfile_potential)
                        
            # indicate that the basis file is already processed
            self.thread.basisfile_potential = ""

        # --- check if there is some new data and if yes, plot it ---
        
        if not self.thread.dataqueue_potential.empty():
        
            # --- storage that allows to draw the hole buffer at once, at least if it is no very large ---
            x = [np.array([])]*3
            y = [np.array([])]*3
            
            while not self.thread.dataqueue_potential.empty() and dataamount < 5000:
                # --- load eigenvalues (energies, y value) and eigenvectors (basis) ---
                filestep, blocknumber, filename = self.thread.dataqueue_potential.get()
                eigensystem = Eigensystem(filename)
                energies = eigensystem.energies
                basis = eigensystem.basis
                
                symmetry = blocknumber % 3
            
                # --- determine which basis elements are within the energy range ---
                boolarr = np.ones(len(energies),dtype=np.bool)
                boolarr[np.isnan(energies)] = False
                energies[np.isnan(energies)] = 0
                
                if self.minE_potential is not None: boolarr &= energies >= self.minE_potential
                if self.maxE_potential is not None: boolarr &= energies <= self.maxE_potential
            
                # cut the energies
                energies = energies[boolarr]
            
                # convert the energies
                energies *= self.converter_y
            
                # cut the basis
                idxarr, = np.nonzero(boolarr)
                transformator = sparse.coo_matrix((np.ones_like(idxarr),(idxarr,np.arange(len(idxarr)))),shape=(basis.shape[1], len(idxarr))).tocsr()
                basis *= transformator
                
                # probability to be in a certain state
                probs = np.abs(basis).power(2) # nState, nBasis
            
                # --- calculate the position (x value) ---
                if self.constDistance and not self.constEField:
                    position = self.get1DPosition(float(eigensystem.params["Ex"]),float(eigensystem.params["Ey"]),float(eigensystem.params["Ez"]))*self.converter_x_potential
                elif self.constDistance and not self.constBField:
                    position = self.get1DPosition(float(eigensystem.params["Bx"]),float(eigensystem.params["By"]),float(eigensystem.params["Bz"]))*self.converter_x_potential
                else:
                    position = float(eigensystem.params["R"])*self.converter_x_potential
            
                # --- draw labels at the beginning of the plotting --- 
                if self.ui.groupbox_plot_labels.isChecked() and filestep == 0:

                    # probability to find a label inside a basis element
                    if self.labelprob_potential is None:
                        self.labelprob_potential = (probs.T*self.labelmat_potential).tocsr() #nBasis, nLabels
                        self.labelprob_energy_potential = [energies]
                        self.labelprob_num_potential = 1
                    else:
                        csr_vappend(self.labelprob_potential,(probs.T*self.labelmat_potential).tocsr())
                        self.labelprob_energy_potential.append(energies)
                        self.labelprob_num_potential += 1
                                                               
                    if not self.labelprob_num_potential < self.thread.numBlocks_potential:
                    
                        self.labelprob_energy_potential = np.concatenate(self.labelprob_energy_potential )
                    
                        # total probability to find a label
                        cumprob = np.array(self.labelprob_potential.sum(axis=0).flat)
                        boolarr = cumprob > 0.1
                                                     
                        # normalize in such a way that the total probability is always one
                        idxarr, = np.nonzero(boolarr)
                        normalizer = sparse.coo_matrix((1/cumprob[idxarr],(idxarr,np.arange(len(idxarr)))),shape=(self.labelprob_potential.shape[1], len(idxarr))).tocsr()
                
                        # store the energetic expectation value of the labels
                        labelenergies = ((self.labelprob_potential*normalizer).T*self.labelprob_energy_potential)
                
                        # store the position of the labels
                        labelposition = position

                        # get size and alpha value of labels
                        size = '{}'.format(max(int(round(self.ui.spinbox_plot_szLabel.value()*11)),1))
                        alpha = int(round(self.ui.spinbox_plot_transpLabel.value()*255))
                
                        # draw the labels
                        for labelstate, labelenergy in zip(self.labelstates_potential[boolarr],labelenergies):
                            sn1, sl1, sj1, sn2, sl2, sj2 = labelstate
                        
                            if self.leftSmallerRight:
                                anchorX = 0
                            else:
                                anchorX = 1

                            text = pg.TextItem(html='<div style="text-align: center; font-size: '+size+'pt;"><span style="color: rgba(0,0,0,255);">'\
                                +'{}{}<sub style="font-size: '.format(int(sn1),self.momentumstrings_potential[int(sl1)]) \
                                +size+'pt;">{}/2</sub>'.format(int(2*sj1))\
                                +' {}{}<sub style="font-size: '.format(int(sn2),self.momentumstrings_potential[int(sl2)]) \
                                +size+'pt;">{}/2</sub>'.format(int(2*sj2))\
                                +'</span></div>',anchor=(anchorX, 0.5),fill=(250,235,215,alpha),border=(255,228,181,255))
                            text.setPos(labelposition, labelenergy)
                            text.setZValue(15)
                            self.ui.graphicsview_potential_plot.addItem(text)
                        
                        posx = labelposition*np.ones_like(labelenergies)+1e-12
                        posy = labelenergies+1e-12
                        curve = PointsItem(np.append(posx, posx-2e-12), np.append(posy, posy-2e-12), 0, 0, (255,255,255))
                        curve.setZValue(5)
                        self.ui.graphicsview_potential_plot.addItem(curve)
                
                        # drawing labels take some time
                        dataamount += 3000
                
                # --- draw colormap ---
                if self.ui.groupbox_plot_overlap.isChecked():
                    # get size and alpha value
                    size = self.ui.spinbox_plot_szOverlap.value()
                    alpha = self.ui.spinbox_plot_transpOverlap.value()
                    
                    # get resolution
                    res = self.ui.spinbox_plot_resolution.value()
                    
                    """if blocknumber not in self.buffer_energiesMap_potential.keys():
                        self.buffer_energiesMap_potential[blocknumber] = {}
                        self.buffer_positionsMap_potential[blocknumber] = {}
                        self.buffer_overlapMap_potential[blocknumber] = {}"""
                    
                    if filestep not in self.buffer_positionsMap_potential.keys():
                        self.buffer_positionsMap_potential[filestep] = position # TODO aufpassen mit submatrix idx
                        self.buffer_energiesMap_potential[filestep] = []
                        self.buffer_overlapMap_potential[filestep] = []
                
                    self.buffer_energiesMap_potential[filestep].append(energies) # TODO aufpassen mit submatrix idx
                    
                    if self.stateidx_potential < 0:
                        self.buffer_overlapMap_potential[filestep].append(np.zeros_like(energies))
                    else:
                        self.buffer_overlapMap_potential[filestep].append(probs[self.stateidx_potential].toarray().flatten()) # TODO aufpassen mit submatrix idx
                                                                                
                    while self.colormap_buffer_minIdx in self.buffer_positionsMap_potential.keys() and self.colormap_buffer_minIdx+1 in self.buffer_positionsMap_potential.keys() and self.colormap_buffer_minIdx+2 in self.buffer_positionsMap_potential.keys():
                        if len(self.buffer_energiesMap_potential[self.colormap_buffer_minIdx]) < self.thread.numBlocks_potential: break
                        if len(self.buffer_energiesMap_potential[self.colormap_buffer_minIdx+1]) < self.thread.numBlocks_potential: break
                        if len(self.buffer_energiesMap_potential[self.colormap_buffer_minIdx+2]) < self.thread.numBlocks_potential: break
                        
                        
                        
                        pos1 = self.buffer_positionsMap_potential[self.colormap_buffer_minIdx]
                        pos2 = self.buffer_positionsMap_potential[self.colormap_buffer_minIdx+1]
                        pos3 = self.buffer_positionsMap_potential[self.colormap_buffer_minIdx+2]
                        energy = np.concatenate(self.buffer_energiesMap_potential[self.colormap_buffer_minIdx+1])
                        overlap = np.concatenate(self.buffer_overlapMap_potential[self.colormap_buffer_minIdx+1])
                        
                        
                        
                        if len(energy) == 0:
                            del self.buffer_energiesMap_potential[self.colormap_buffer_minIdx]
                            del self.buffer_positionsMap_potential[self.colormap_buffer_minIdx]
                            del self.buffer_overlapMap_potential[self.colormap_buffer_minIdx]
                            self.colormap_buffer_minIdx += 1
                                
                            continue
                        
                        
                                                
                        #posMin = (pos1+pos2)/2 # TODO schauen wie es sich mit der Plotreihenfolge vertraegt
                        #posMax = (pos2+pos3)/2
                        
                        posMin = pos1
                        posMax = pos3
                        
                        
                        

                        
                        
                        
                        
                        
                        
                        # first call
                        if np.isnan(self.displayunits2pixelunits_x) or np.isnan(self.displayunits2pixelunits_y):
                            
                            if self.maxE_potential is not None:
                                energyMax = self.maxE_potential*self.converter_y
                            else:
                                energyMax = np.nanmax(energy)
                                
                            if self.minE_potential is not None:
                                energyMin = self.minE_potential*self.converter_y
                            else:
                                energyMin = np.nanmin(energy)
                        
                            if posMax == posMin or energyMax == energyMin :
                                del self.buffer_energiesMap_potential[self.colormap_buffer_minIdx]
                                del self.buffer_positionsMap_potential[self.colormap_buffer_minIdx]
                                del self.buffer_overlapMap_potential[self.colormap_buffer_minIdx]
                                self.colormap_buffer_minIdx += 1
                                
                                continue
                            
                            # TODO get resolution of plot area
                            
                            
                            
                            self.height_pixelunits = res
                            
                            enlargement = max(np.round((self.height_pixelunits/self.steps-2)/2),0)
                            
                            self.width_pixelunits = 5+4*enlargement #np.round((posMax-posMin)*self.displayunits2pixelunits_x)
                            
                            
                            self.displayunits2pixelunits_x = self.width_pixelunits/(posMax-posMin)
                            self.displayunits2pixelunits_y = self.height_pixelunits/(energyMax-energyMin)
                            self.smootherX = (enlargement*2+2)*1/2
                            self.smootherY = self.height_pixelunits*1/150*size #250/self.steps
                            
                        
                            #margin_pixelunits = np.round(max(self.width_pixelunits/2,3*self.smootherY))
                            #margin_displayunits = margin_pixelunits/self.displayunits2pixelunits_y
                        
                            
                            #height_displayunits = self.height_pixelunits/self.displayunits2pixelunits_y
                            
                            
                            
                            
                            """self.displayunits2pixelunits_x = (5+4)/(posMax-posMin)
                            self.displayunits2pixelunits_y = (5+4)*self.steps/(energyMax-energyMin)
                            self.smootherX = 2
                            self.smootherY = self.steps/17
                            
                            
                        
                            width_pixelunits = 9 #np.round((posMax-posMin)*self.displayunits2pixelunits_x)
                        
                            margin_pixelunits = np.round(max(width_pixelunits/2,3*self.smootherY))
                            margin_displayunits = margin_pixelunits/self.displayunits2pixelunits_y
                        
                            self.height_pixelunits = np.round((energyMax-energyMin)*self.displayunits2pixelunits_y)+2*margin_pixelunits
                            height_displayunits = self.height_pixelunits/self.displayunits2pixelunits_y"""
                        
                            self.yMin = energyMin#-margin_displayunits
                            self.yMax = energyMax#self.yMin+height_displayunits
                        
                        #width_pixelunits = np.round((posMax-posMin)*self.displayunits2pixelunits_x)
                                                
                        map = np.zeros((self.width_pixelunits,self.height_pixelunits)) # x-y-coordinate system, origin is at the bottom left corner 
                        
                                             
                        
                        """map[::,::2] = 3
                        map[-1,-1] = 10
                        map[0,0] = 10"""
                        
                        posLeft = (pos1+pos2)/2
                        posRight = (pos2+pos3)/2
                        
                        idx_position = bytescale(np.array([pos1, pos2, pos3]), low=0, high=self.width_pixelunits-1, cmin=pos1, cmax=pos3)
                        idx_left = bytescale(posLeft, low=0, high=self.width_pixelunits-1, cmin=pos1, cmax=pos3)
                        idx_right = bytescale(posRight, low=0, high=self.width_pixelunits-1, cmin=pos1, cmax=pos3)
                        
                        
                        for i in range(3):
                            energy = np.concatenate(self.buffer_energiesMap_potential[self.colormap_buffer_minIdx+i])
                            overlap = np.concatenate(self.buffer_overlapMap_potential[self.colormap_buffer_minIdx+i])
                            
                            boolarr = (energy >= self.yMin) & (energy <= self.yMax)
                            energy = energy[boolarr]
                            overlap = overlap[boolarr]
                            
                            if self.logscale:
                                overlap[overlap < 1e-2] = 1e-2
                                overlap = (2+np.log10(overlap))/2
                            
                            idx_energy = bytescale(energy, low=0, high=self.height_pixelunits-1, cmin=self.yMin, cmax=self.yMax)
                            mapmat = sparse.coo_matrix((overlap,(idx_energy,np.arange(len(idx_energy)))),shape=(self.height_pixelunits, len(idx_energy))).sum(axis=-1)
                            
                            map[idx_position[i],:] = mapmat.flat

                            
                            
                        
                        
                        """print(idx_position, idx_left, idx_right)
                        
                        # 0, 3, 6
                        
                        
                        idx_energy = bytescale(energy, low=0, high=height_pixelunits-1, cmin=yMin, cmax=yMax)
                        mapmat = sparse.coo_matrix((overlap,(idx_energy,np.arange(len(idx_energy)))),shape=(height_pixelunits, len(idx_energy))).sum(axis=-1)
                        
                        
                        
                        map[(width_pixelunits-1)/2,:] = mapmat.flat"""
                        
                        map = gaussian_filter(map,(self.smootherX,self.smootherY),mode='mirror')
                        
                        map = map[idx_left:idx_right]
                        
                        normalizer = np.zeros((self.width_pixelunits,2*3*self.smootherY+1))
                        normalizer[(self.width_pixelunits-1)/2,3*self.smootherY] = 1
                        normalizer = gaussian_filter(normalizer,(self.smootherX,self.smootherY),mode='mirror')
                        
                        map /= np.max(normalizer)
                        
                        # TODO no auto limits
                        
                        
                        
                        #height = int(round(width*self.ratio_pixelunits_potential/self.ratio_displayunits_potential*ratio + 2*margin_pixelunits))
                        
                        
                        
                        
                        #ratio = (energyMax-energyMin)/(posMax-posMin)
                        
                       
                        
                        
                        
                        #height = int(round(width*self.ratio_pixelunits_potential/self.ratio_displayunits_potential*ratio + 2*margin_pixelunits))
                                                
                        
                        #0*marginInPointunit*(energyMax-energyMin)/map.shape[1]
                        
                        """map = np.zeros((width_pixelunits,height_pixelunits)) # x-y-coordinate system, origin is at the bottom left corner                        
                        
                        map[::,::2] = 3
                        map[-1,-1] = 10
                        map[0,0] = 10
                        
                        
                        
                        aspectRatio = map.shape[1]/map.shape[0] # TODO calculate it before"""
                        
                        
                        
                        
                        
                        img = pg.ImageItem(image=map, opacity=alpha, autoDownsample=True, lut=self.lut, levels=[0,1])
                        img.setRect(QtCore.QRectF(posLeft-0.5/self.displayunits2pixelunits_x,self.yMin-0.5/self.displayunits2pixelunits_y,posRight-posLeft,self.yMax-self.yMin+1/self.displayunits2pixelunits_y)) # xMin, yMin, xSize, ySize # TODO energyMin anpassen wegen Pixelgroesse
                        img.setZValue(3)
                        self.ui.graphicsview_potential_plot.addItem(img)
                        
                        dataamount += len(energy)*10
                    
                        del self.buffer_energiesMap_potential[self.colormap_buffer_minIdx]
                        del self.buffer_positionsMap_potential[self.colormap_buffer_minIdx]
                        del self.buffer_overlapMap_potential[self.colormap_buffer_minIdx]
                        self.colormap_buffer_minIdx += 1
                        
                
                    
                    
                    
                                
                # --- draw lines ---
                if self.ui.groupbox_plot_lines.isChecked():
                    if blocknumber not in self.buffer_basis_potential.keys():
                        self.buffer_basis_potential[blocknumber] = {}
                        self.buffer_energies_potential[blocknumber] = {}
                        self.buffer_positions_potential[blocknumber] = {}
                        self.lines_buffer_minIdx[blocknumber] = 0
                
                    self.buffer_basis_potential[blocknumber][filestep] = basis
                    self.buffer_energies_potential[blocknumber][filestep] = energies
                    self.buffer_positions_potential[blocknumber][filestep] = position
                                            
                    # get size and alpha value of points
                    size = self.ui.spinbox_plot_szLine.value()
                    alpha = self.ui.spinbox_plot_transpLine.value()*255
                                        
                    while self.lines_buffer_minIdx[blocknumber] in self.buffer_basis_potential[blocknumber].keys() and self.lines_buffer_minIdx[blocknumber]+1 in self.buffer_basis_potential[blocknumber].keys():
                        # determine the data to plot
                        overlap = np.abs(self.buffer_basis_potential[blocknumber][self.lines_buffer_minIdx[blocknumber]].T*self.buffer_basis_potential[blocknumber][self.lines_buffer_minIdx[blocknumber]+1]) # nBasis first, nBasis second
                        overlap.data[overlap.data <= np.sqrt(0.5)] = 0
                        overlap.eliminate_zeros()
                        overlap = overlap.tocoo()
                            
                        iFirst = overlap.row
                        iSecond = overlap.col
                        
                        ydata = np.transpose([self.buffer_energies_potential[blocknumber][self.lines_buffer_minIdx[blocknumber]][iFirst],self.buffer_energies_potential[blocknumber][self.lines_buffer_minIdx[blocknumber]+1][iSecond]])
                        xdata = np.ones_like(ydata)
                        xdata[:,0] *= self.buffer_positions_potential[blocknumber][self.lines_buffer_minIdx[blocknumber]]
                        xdata[:,1] *= self.buffer_positions_potential[blocknumber][self.lines_buffer_minIdx[blocknumber]+1]
                        
                        if len(iFirst) != 0:
                            # plot the data
                            curve = MultiLine(xdata, ydata, size, alpha, self.symmetrycolors[symmetry])
                            curve.setZValue(4)
                            self.ui.graphicsview_potential_plot.addItem(curve)
                        
                        del self.buffer_basis_potential[blocknumber][self.lines_buffer_minIdx[blocknumber]]
                        del self.buffer_energies_potential[blocknumber][self.lines_buffer_minIdx[blocknumber]]
                        del self.buffer_positions_potential[blocknumber][self.lines_buffer_minIdx[blocknumber]]
                        
                        self.lines_buffer_minIdx[blocknumber] += 1
                        dataamount += len(iFirst)*10
                    
                # --- store data to plot several points at once ---
                if self.ui.groupbox_plot_points.isChecked():
                    x[symmetry] = np.append(x[symmetry],position*np.ones_like(energies))
                    y[symmetry] = np.append(y[symmetry],energies)
                    
                    dataamount += len(x)            
        
            # --- plot the stored data ---
            if self.ui.groupbox_plot_points.isChecked() and len(x) > 0:
        
                # get size and alpha value of points
                size = self.ui.spinbox_plot_szPoint.value()
                alpha = self.ui.spinbox_plot_transpPoint.value()*255
                
                # plot the basis elements
                
                # loop over symmetries
                for s in range(3):
                    if len(x[s]) == 0: continue
                    curve = PointsItem(x[s], y[s], size, alpha, self.symmetrycolors[s])
                    curve.setZValue(10)
                    self.ui.graphicsview_potential_plot.addItem(curve)
        
            # --- update the graphics view ---
            self.ui.graphicsview_potential_plot.repaint()
        
        # === process field maps ===
             
        for idx in range(2):
        
            dataqueue = [self.thread.dataqueue_field1, self.thread.dataqueue_field2][idx]
            basisfile = [self.thread.basisfile_field1, self.thread.basisfile_field2][idx]
            graphicsview_plot = [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot]
            
            if not self.thread.samebasis:
                minE = [self.minE_field1, self.minE_field2][idx]
                maxE = [self.maxE_field1, self.maxE_field2][idx]
            else:
                if self.minE_field1 is None or self.minE_field2 is None: minE = None
                else: minE = min(self.minE_field1, self.minE_field2)
                if self.maxE_field1 is None or self.maxE_field2 is None: maxE = None
                else: maxE = max(self.maxE_field1, self.maxE_field2)
            
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
                self.labelmat[idx] = sparse.coo_matrix((np.ones_like(cumsum),(np.arange(len(cumsum)),cumsum-1)),shape=(len(cumsum), len(self.labelstates[idx]))).tocsr() #nStates, nLabels
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
            
                while not dataqueue.empty() and dataamount < 5000: # stop loop if enough data is collected
                
                    # --- load eigenvalues (energies, y value) and eigenvectors (basis) ---
                    filestep, blocknumber, filename = dataqueue.get()
                    eigensystem = Eigensystem(filename)
                    energies = eigensystem.energies
                    basis = eigensystem.basis
                
                    # --- determine which basis elements are within the energy range ---
                    boolarr = np.ones(len(energies),dtype=np.bool)
                    boolarr[np.isnan(energies)] = False
                    energies[np.isnan(energies)] = 0
                
                    if minE is not None: boolarr &= energies >= minE
                    if maxE is not None: boolarr &= energies <= maxE
                
                    # cut the energies
                    energies = energies[boolarr]
                
                    # convert the energies
                    energies *= self.converter_y
                
                    # cut the basis
                    idxarr, = np.nonzero(boolarr)
                    transformator = sparse.coo_matrix((np.ones_like(idxarr),(idxarr,np.arange(len(idxarr)))),shape=(basis.shape[1], len(idxarr))).tocsr()
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
                        position = self.get1DPosition(float(eigensystem.params["Bx"]),float(eigensystem.params["By"]),float(eigensystem.params["Bz"]))*self.converter_x_field
                    else:
                        position = self.get1DPosition(float(eigensystem.params["Ex"]),float(eigensystem.params["Ey"]),float(eigensystem.params["Ez"]))*self.converter_x_field

                    # --- draw labels at the beginning of the plotting --- 
                    if self.ui.groupbox_plot_labels.isChecked() and filestep == 0:
                
                        # probability to find a label inside a basis element
                        labelprob = probs.T*self.labelmat[idx] #nBasis, nLabels
                    
                        # total probability to find a label
                        cumprob = np.array(labelprob.sum(axis=0).flat)
                        boolarr = cumprob > 0.1
                    
                        # normalize in such a way that the total probability is always one
                        idxarr, = np.nonzero(boolarr)
                        normalizer = sparse.coo_matrix((1/cumprob[idxarr],(idxarr,np.arange(len(idxarr)))),shape=(labelprob.shape[1], len(idxarr))).tocsr()
                    
                        # store the energetic expectation value of the labels
                        labelenergies = ((labelprob*normalizer).T*energies)
                        
                        if len(labelenergies) == 0: continue
                    
                        # store the position of the labels
                        labelposition = position

                        # get size and alpha value of labels
                        size = '{}'.format(max(int(round(self.ui.spinbox_plot_szLabel.value()*11)),1))
                        alpha = int(round(self.ui.spinbox_plot_transpLabel.value()*255))
                    
                        # draw the labels
                        for labelstate, labelenergy in zip(self.labelstates[idx][boolarr],labelenergies):
                            sn, sl, sj = labelstate
                            
                            if self.leftSmallerRight:
                                anchorX = 0
                            else:
                                anchorX = 1
                            
                            for j in range(2):
                                if not self.thread.samebasis and j != idx: continue
                                text = pg.TextItem(html='<div style="text-align: center; font-size: '+size+'pt;">'\
                                    +'<span style="color: rgba(0,0,0,255);">{}{}<sub style="font-size: '.format(int(sn),self.momentumstrings[idx][int(sl)]) \
                                    +size+'pt;">{}/2</sub></span></div>'.format(int(2*sj)),anchor=(anchorX, 0.5),fill=(250,235,215,alpha),border=(255,228,181,255))                    
                                text.setPos(labelposition, labelenergy)
                                text.setZValue(15)
                                graphicsview_plot[j].addItem(text)
                            
                        for j in range(2):
                            if not self.thread.samebasis and j != idx: continue
                            posx = labelposition*np.ones_like(labelenergies)+1e-12
                            posy = labelenergies+1e-12
                            curve = PointsItem(np.append(posx, posx-2e-12), np.append(posy, posy-2e-12), 0, 0, (255,255,255))
                            curve.setZValue(5)
                            graphicsview_plot[j].addItem(curve)               
                    
                        # drawing labels take some time
                        dataamount += 3000
                    
                    # --- draw lines ---
                    if self.ui.groupbox_plot_lines.isChecked():
                        self.buffer_basis[idx][filestep] = basis
                        self.buffer_energies[idx][filestep] = energies
                        self.buffer_positions[idx][filestep] = position
                        self.buffer_boolarr[idx][filestep] = []
                        
                        # cut the momenta to reasonable values
                        momentum[momentum>len(self.momentumcolors)-2] = len(self.momentumcolors)-2
                        momentum[momentum<0] = len(self.momentumcolors)-1
                        
                        # loop over momenta
                        for i in range(len(self.momentumcolors)):
                            # determine which basis elements have the current momentum
                            boolarr = momentum == i
                            self.buffer_boolarr[idx][filestep].append(boolarr)
                                                
                        # get size and alpha value of points
                        size = self.ui.spinbox_plot_szLine.value()
                        alpha = self.ui.spinbox_plot_transpLine.value()*255
                                                
                        while self.lines_buffer_minIdx_field in self.buffer_basis[idx].keys() and self.lines_buffer_minIdx_field+1 in self.buffer_basis[idx].keys():
                            # determine the data to plot
                            overlap = np.abs(self.buffer_basis[idx][self.lines_buffer_minIdx_field].T*self.buffer_basis[idx][self.lines_buffer_minIdx_field+1]) # nBasis first, nBasis second
                            overlap.data[overlap.data <= np.sqrt(0.8)] = 0
                            overlap.eliminate_zeros()
                            overlap = overlap.tocoo()
                                
                            iFirst = overlap.row
                            iSecond = overlap.col
                            
                            ydata = np.transpose([self.buffer_energies[idx][self.lines_buffer_minIdx_field][iFirst],self.buffer_energies[idx][self.lines_buffer_minIdx_field+1][iSecond]])
                            xdata = np.ones_like(ydata)
                            xdata[:,0] *= self.buffer_positions[idx][self.lines_buffer_minIdx_field]
                            xdata[:,1] *= self.buffer_positions[idx][self.lines_buffer_minIdx_field+1]
                        
                            # loop over momenta
                            for i in range(len(self.momentumcolors)):
                                if np.sum(self.buffer_boolarr[idx][self.lines_buffer_minIdx_field][i]) == 0: continue
                                boolarr = self.buffer_boolarr[idx][self.lines_buffer_minIdx_field][i][iFirst]
                                if np.sum(boolarr) == 0: continue
                            
                                # determine the associated color
                                color = self.momentumcolors[i]
                                
                                # plot the data
                                for j in range(2):
                                    if not self.thread.samebasis and j != idx: continue
                                    curve = MultiLine(xdata[boolarr], ydata[boolarr], size, alpha, color) # TODO alpha and color der Funktion zusammen uebergeben
                                    curve.setZValue(4)
                                    graphicsview_plot[j].addItem(curve)
                            
                            del self.buffer_basis[idx][self.lines_buffer_minIdx_field]
                            del self.buffer_energies[idx][self.lines_buffer_minIdx_field]
                            del self.buffer_positions[idx][self.lines_buffer_minIdx_field]
                            del self.buffer_boolarr[idx][self.lines_buffer_minIdx_field]
                            
                            self.lines_buffer_minIdx_field += 1
                            
                            dataamount += len(iFirst)*10
                            
                            
                            
                    
                    # --- store data to plot several points at once ---
                    if self.ui.groupbox_plot_points.isChecked():
                        x = np.append(x,position*np.ones_like(energies))
                        y = np.append(y,energies)
                        l = np.append(l,momentum)
                        
                        dataamount += len(x)
                    
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
            
            # Delete buffers
            self.buffer_basis = [{}]*2
            self.buffer_energies = [{}]*2
            self.buffer_positions = [{}]*2
            self.buffer_boolarr = [{}]*2
            self.buffer_basis_potential = {}
            self.buffer_energies_potential = {}
            self.buffer_positions_potential = {}
            
            self.buffer_energiesMap = [{}]*2
            self.buffer_positionsMap = [{}]*2
            self.buffer_energiesMap_potential = {}
            self.buffer_positionsMap_potential = {}
            self.buffer_overlapMap_potential = {}
            
            self.labelprob_potential = None
            self.labelprob_num_potential = 0
            
            self.lines_buffer_minIdx = {}
            self.colormap_buffer_minIdx = 0
            self.lines_buffer_minIdx_field = 0
        
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
    
    @QtCore.pyqtSlot(bool) # TODO
    def toggleAntialiasing(self):
        checked = self.ui.checkbox_plot_antialiasing.isChecked()
        for plotarea in [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot, self.ui.graphicsview_potential_plot]:
            plotarea.setAntialiasing(checked)
        pg.setConfigOptions(antialias=checked) # TODO
    
    @QtCore.pyqtSlot(bool) # TODO
    def togglePairbasis(self):
        checked = self.ui.radiobutton_system_pairbasisDefined.isChecked()
        self.ui.widget_system_pair.setEnabled(checked)
    
    @QtCore.pyqtSlot(bool) # TODO
    def toggleOverlapstate(self):
        checked = self.ui.radiobutton_plot_overlapDefined.isChecked()
        self.ui.widget_plot_qn.setEnabled(checked)
    
    @QtCore.pyqtSlot(bool) # TODO
    def toggleYScale(self):
        log = self.ui.radiobutton_plot_log.isChecked()
        if log:
            self.ui.label_plot_1st.setText("< 0.01") # TODO
            self.ui.label_plot_2nd.setText("0.1") # TODO
            self.logscale = True
        else:
            self.ui.label_plot_1st.setText("0")
            self.ui.label_plot_2nd.setText("0.5")
            self.logscale = False
    
    @QtCore.pyqtSlot(str) # TODO
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
                    params = self.systemdict.paramsInOriginalunits()
                    json.dump(params, f, indent=4, sort_keys=True)
                    
                with open(self.path_plot_last, 'w') as f:
                    params = self.plotdict.paramsInOriginalunits()
                    params["gradientwidget"] = self.ui.gradientwidget_plot_gradient.saveState()
                    json.dump(params, f, indent=4, sort_keys=True)
                    
                with open(self.path_tabs_last, 'w') as f:
                    params = dict()
                    params["config"] = self.ui.tabwidget_config.currentIndex()
                    params["plotter"] = self.ui.tabwidget_plotter.currentIndex()
                    params["system"] = self.ui.toolbox_system.currentIndex()
                    json.dump(params, f, indent=4, sort_keys=True)
            
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
                #self.sameSpecies = self.getSameSpecies()
            
                if self.senderbutton in [self.ui.pushbutton_field1_calc, self.ui.pushbutton_potential_calc] or self.samebasis:
                    self.ui.graphicsview_field1_plot.clear()
                    self.ui.graphicsview_field1_plot.setLabel('left', 'Energy ('+str(Units.energy)+')')
                    if self.constEField and not self.constBField:
                        self.ui.graphicsview_field1_plot.setLabel('bottom', 'Magnetic field ('+str(Units.bfield)+')')
                        self.converter_x_field = Converter.fromAU(1,Units.bfield).magnitude
                        posMin = self.get1DPosition(self.systemdict['minBx'].magnitude,self.systemdict['minBy'].magnitude,self.systemdict['minBz'].magnitude)
                        posMax = self.get1DPosition(self.systemdict['maxBx'].magnitude,self.systemdict['maxBy'].magnitude,self.systemdict['maxBz'].magnitude)
                    else:
                        self.ui.graphicsview_field1_plot.setLabel('bottom', 'Electric field ('+str(Units.efield)+')')
                        self.converter_x_field = Converter.fromAU(1,Units.efield).magnitude
                        posMin = self.get1DPosition(self.systemdict['minEx'].magnitude,self.systemdict['minEy'].magnitude,self.systemdict['minEz'].magnitude)
                        posMax = self.get1DPosition(self.systemdict['maxEx'].magnitude,self.systemdict['maxEy'].magnitude,self.systemdict['maxEz'].magnitude)
                if self.senderbutton in [self.ui.pushbutton_field2_calc, self.ui.pushbutton_potential_calc] or self.samebasis:
                    self.ui.graphicsview_field2_plot.clear()
                    self.ui.graphicsview_field2_plot.setLabel('left', 'Energy ('+str(Units.energy)+')')
                    if self.constEField and not self.constBField:
                        self.ui.graphicsview_field2_plot.setLabel('bottom', 'Magnetic field ('+str(Units.bfield)+')')
                        self.converter_x_field = Converter.fromAU(1,Units.bfield).magnitude
                        posMin = self.get1DPosition(self.systemdict['minBx'].magnitude,self.systemdict['minBy'].magnitude,self.systemdict['minBz'].magnitude)
                        posMax = self.get1DPosition(self.systemdict['maxBx'].magnitude,self.systemdict['maxBy'].magnitude,self.systemdict['maxBz'].magnitude)
                    else:
                        self.ui.graphicsview_field2_plot.setLabel('bottom', 'Electric field ('+str(Units.efield)+')')
                        self.converter_x_field = Converter.fromAU(1,Units.efield).magnitude
                        posMin = self.get1DPosition(self.systemdict['minEx'].magnitude,self.systemdict['minEy'].magnitude,self.systemdict['minEz'].magnitude)
                        posMax = self.get1DPosition(self.systemdict['maxEx'].magnitude,self.systemdict['maxEy'].magnitude,self.systemdict['maxEz'].magnitude)
                if self.senderbutton == self.ui.pushbutton_potential_calc:
                    self.ui.graphicsview_potential_plot.clear()
                    self.ui.graphicsview_potential_plot.setLabel('left', 'Energy ('+str(Units.energy)+')')
                    if self.constDistance and not self.constEField:
                        self.ui.graphicsview_potential_plot.setLabel('bottom', 'Electric field ('+str(Units.efield)+')')
                        self.converter_x_potential = Converter.fromAU(1,Units.efield).magnitude
                        posMin = self.get1DPosition(self.systemdict['minEx'].magnitude,self.systemdict['minEy'].magnitude,self.systemdict['minEz'].magnitude)
                        posMax = self.get1DPosition(self.systemdict['maxEx'].magnitude,self.systemdict['maxEy'].magnitude,self.systemdict['maxEz'].magnitude)
                    elif self.constDistance and not self.constBField:
                        self.ui.graphicsview_potential_plot.setLabel('bottom', 'Magnetic field ('+str(Units.bfield)+')')
                        self.converter_x_potential = Converter.fromAU(1,Units.bfield).magnitude
                        posMin = self.get1DPosition(self.systemdict['minBx'].magnitude,self.systemdict['minBy'].magnitude,self.systemdict['minBz'].magnitude)
                        posMax = self.get1DPosition(self.systemdict['maxBx'].magnitude,self.systemdict['maxBy'].magnitude,self.systemdict['maxBz'].magnitude)
                    else:
                        self.ui.graphicsview_potential_plot.setLabel('bottom', 'Interatomic distance ('+str(Units.length)+')')
                        self.converter_x_potential = Converter.fromAU(1,Units.length).magnitude
                        posMin = self.systemdict['minR'].magnitude
                        posMax = self.systemdict['maxR'].magnitude
                        
                self.converter_y = Converter.fromAU(1,Units.energy).magnitude
                self.leftSmallerRight = posMin < posMax
                self.steps = self.systemdict['steps']
                
                if self.ui.groupbox_plot_overlap.isChecked():
                    if self.plotdict["overlapUnperturbed"]:
                        n1 = self.systemdict['n1']
                        l1 = self.systemdict['l1']
                        j1 = self.systemdict['j1']
                        m1 = self.systemdict['m1']
                        n2 = self.systemdict['n2']
                        l2 = self.systemdict['l2']
                        j2 = self.systemdict['j2']
                        m2 = self.systemdict['m2']
                    else:
                        n1 = self.plotdict['n1']
                        l1 = self.plotdict['l1']
                        j1 = self.plotdict['j1']
                        m1 = self.plotdict['m1']
                        n2 = self.plotdict['n2']
                        l2 = self.plotdict['l2']
                        j2 = self.plotdict['j2']
                        m2 = self.plotdict['m2']
                    self.overlapstate = np.array([n1,l1,j1,m1,n2,l2,j2,m2])
                    
                self.lut = self.ui.gradientwidget_plot_gradient.getLookupTable(512)
                self.lut[0] = [255,255,255]
                    
                # Set limits
                self.minE_field1 = self.plotdict["minE_field1"]
                if self.minE_field1 is not None:
                    self.minE_field1 = Converter.toAU(self.minE_field1).magnitude
                    self.ui.graphicsview_field1_plot.setLimits(yMin = self.minE_field1 * self.converter_y)
                else:
                    self.ui.graphicsview_field1_plot.setLimits(yMin = None)
                    
                self.minE_field2 = self.plotdict["minE_field2"]
                if self.minE_field2 is not None:
                    self.minE_field2 = Converter.toAU(self.minE_field2).magnitude
                    self.ui.graphicsview_field2_plot.setLimits(yMin = self.minE_field2 * self.converter_y)
                else:
                    self.ui.graphicsview_field2_plot.setLimits(yMin = None)
                
                self.minE_potential = self.plotdict["minE_potential"]
                if self.minE_potential is not None:
                    self.minE_potential = Converter.toAU(self.minE_potential).magnitude
                    self.ui.graphicsview_potential_plot.setLimits(yMin = self.minE_potential * self.converter_y)
                else:
                    self.ui.graphicsview_potential_plot.setLimits(yMin = None)
                
                self.maxE_field1 = self.plotdict["maxE_field1"]
                if self.maxE_field1 is not None:
                    self.maxE_field1 = Converter.toAU(self.maxE_field1).magnitude
                    self.ui.graphicsview_field1_plot.setLimits(yMax = self.maxE_field1 * self.converter_y)
                else:
                    self.ui.graphicsview_field1_plot.setLimits(yMax = None)
                    
                self.maxE_field2 = self.plotdict["maxE_field2"]
                if self.maxE_field2 is not None:
                    self.maxE_field2 = Converter.toAU(self.maxE_field2).magnitude
                    self.ui.graphicsview_field2_plot.setLimits(yMax = self.maxE_field2 * self.converter_y)
                else:
                    self.ui.graphicsview_field2_plot.setLimits(yMax = None)
                
                self.maxE_potential = self.plotdict["maxE_potential"]
                if self.maxE_potential is not None:
                    self.maxE_potential = Converter.toAU(self.maxE_potential).magnitude
                    self.ui.graphicsview_potential_plot.setLimits(yMax = self.maxE_potential * self.converter_y)
                else:
                    self.ui.graphicsview_potential_plot.setLimits(yMax = None)
                
                
            
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
                params = self.systemdict.paramsInOriginalunits()
                json.dump(params, f, indent=4, sort_keys=True)
            self.systemfile = filename
            self.systempath = os.path.dirname(filename)
    
    @QtCore.pyqtSlot()
    def savePlotConf(self):
        path = self.plotfile if self.plotfile is not None else self.plotpath
        filename = QtGui.QFileDialog.getSaveFileName(self, \
            "Save plot configuration",path, "json (*.json)")
        
        if filename:
            with open(filename, 'w') as f:
                params = self.plotdict.paramsInOriginalunits()
                params["gradientwidget"] = self.ui.gradientwidget_plot_gradient.saveState()
                json.dump(params, f, indent=4, sort_keys=True)
            self.plotfile = filename
            self.plotpath = os.path.dirname(filename)
    
    @QtCore.pyqtSlot()
    def openSystemConf(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, \
            "Open system configuration",self.systempath, "json (*.json)")
        
        if filename:
            with open(filename, 'r') as f:
                params = json.load(f)
                self.systemdict.load(params)
            self.systemfile = filename
            self.systempath = os.path.dirname(filename)
    
    @QtCore.pyqtSlot()
    def openPlotConf(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, \
            "Open plot configuration",self.systempath, "json (*.json)")
        
        if not (filename == ""):
            with open(filename, 'r') as f:
                params = json.load(f)
                self.ui.gradientwidget_plot_gradient.restoreState(params["gradientwidget"])
                del params["gradientwidget"]
                self.plotdict.load(params)
            self.plotfile = filename
            self.plotpath = os.path.dirname(filename)
            
    def closeEvent(self, event):
        # Kill c++ program if necessary 
        self.abortCalculation()
        
        # Save last settings
        with open(self.path_system_last, 'w') as f:
            params = self.systemdict.paramsInOriginalunits()
            json.dump(params, f, indent=4, sort_keys=True)
        
        with open(self.path_plot_last, 'w') as f:
            params = self.plotdict.paramsInOriginalunits()
            params["gradientwidget"] = self.ui.gradientwidget_plot_gradient.saveState()
            json.dump(params, f, indent=4, sort_keys=True)
        
        with open(self.path_tabs_last, 'w') as f:
            params = dict()
            params["config"] = self.ui.tabwidget_config.currentIndex()
            params["plotter"] = self.ui.tabwidget_plotter.currentIndex()
            params["system"] = self.ui.toolbox_system.currentIndex()
            json.dump(params, f, indent=4, sort_keys=True)
        
        # Close everything
        super().closeEvent(event)

def main():
    app = QtGui.QApplication(sys.argv)
    form = MainWindow()
    form.show()
    app.exec_()
    
if __name__ == "__main__":
    main()
