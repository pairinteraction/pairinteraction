#!/usr/bin/env python3

version_settings = 1
version_cache = 1

import sys
from pint import UnitRegistry
from pint.unit import UndefinedUnitError
from PyQt5 import QtCore, QtGui
from plotter import Ui_plotwindow # pyuic4 plotter.ui > plotter.py or py3uic4 plotter.ui > plotter.py
import pyqtgraph as pg
import pyqtgraph.exporters
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
"""from palettable import cubehelix"""
from operator import itemgetter
import zipfile
from scipy import io
from io import StringIO, BytesIO
from shutil import copyfile
import shutil

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

def csc_happend(a,b):
    """ Takes in 2 csc_matrices and appends the second one to the right of the first one."""

    a.data = np.hstack((a.data,b.data))
    a.indices = np.hstack((a.indices,b.indices))
    a.indptr = np.hstack((a.indptr,(b.indptr + a.nnz)[1:]))
    a._shape = (b.shape[0],a.shape[1]+b.shape[1])

#http://stackoverflow.com/questions/15992857/efficient-way-to-get-the-max-of-each-row-for-large-sparse-matrix
def csr_keepmax(a):
    boolarr = np.diff(a.indptr)>0
    ret = np.maximum.reduceat(a.data, a.indptr[:-1][boolarr])
    a.data[a.data != np.repeat(ret,np.diff(a.indptr)[boolarr])] = 0
    a.eliminate_zeros()


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
        store["overlapDefined"] = {'widget': ui.radiobutton_plot_overlapDefined}
        store["overlapTransition"] = {'widget': ui.radiobutton_plot_overlapTransition}
        store["transition"] = {'widget': ui.combobox_plot_overlapTransition}
        store["connectionthreshold"] = {'widget': ui.spinbox_plot_connectionthreshold}
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
        store["autorange"] = {'widget': ui.checkbox_plot_autorange}

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
        locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
        QtCore.QLocale.setDefault(QtCore.QLocale(locale.getlocale()[0]))
        
        super().__init__(parent)
        
        del pg.graphicsItems.GradientEditorItem.Gradients['greyclip']
        #del pg.graphicsItems.GradientEditorItem.Gradients['grey']
        del pg.graphicsItems.GradientEditorItem.Gradients['cyclic']
        del pg.graphicsItems.GradientEditorItem.Gradients['spectrum']
        del pg.graphicsItems.GradientEditorItem.Gradients['bipolar']
        
        
        """color = cubehelix.Cubehelix.make(sat=1.8,n=7,rotation=1.21,start=1.2,reverse=True).colors[::-1]
        color = np.append(color, [[255]]*len(color), axis=1).astype(np.ubyte)"""
        color = [
            [  0,   0,   0, 255],
            [ 12,  67,   0, 255],
            [  0, 123, 118, 255],
            [122, 109, 240, 255],
            [255, 121, 197, 255],
            [247, 204, 159, 255],
            [255, 255, 255, 255]
        ]
        pos = np.linspace(0,1,len(color))
        pg.graphicsItems.GradientEditorItem.Gradients['cubehelix1'] = {'mode':'rgb', 'ticks':[(p, tuple(c)) for c, p in zip(color,pos)]}
        
        
        """color = cubehelix.Cubehelix.make(sat=1.5,n=7,rotation=-1.0,start=0.9,reverse=True).colors[::-1]
        color = np.append(color, [[255]]*len(color), axis=1).astype(np.ubyte)"""
        color = [
            [  0,   0,   0, 255],
            [ 75,  19,  77, 255],
            [ 63,  80, 167, 255],
            [ 44, 164, 156, 255],
            [117, 206, 113, 255],
            [226, 215, 161, 255],
            [255, 255, 255, 255]
        ]
        pos = np.linspace(0,1,len(color))
        pg.graphicsItems.GradientEditorItem.Gradients['cubehelix2'] = {'mode':'rgb', 'ticks':[(p, tuple(c)) for c, p in zip(color,pos)]}
        
        for k, v in pg.graphicsItems.GradientEditorItem.Gradients.items():
            mapticks = [(1-p, c) for p, c in v['ticks']]
            #mappos = [1-p for p, c in v['ticks']]
            #mapticks[np.argmin(mappos)] = (mapticks[np.argmin(mappos)][0],(245,245,245,255))
            pg.graphicsItems.GradientEditorItem.Gradients[k]['ticks'] = mapticks
        
        
        
        self.ui = Ui_plotwindow()
        self.ui.setupUi(self)
        
        self.invalidQuantumnumbers = [False, False, False, False]
        
        self.samebasis = False
        
        self.systemdict = SystemDict(self.ui)
        self.plotdict = PlotDict(self.ui)
        
        if os.name == 'nt': self.userpath = os.path.expanduser('~user')
        else: self.userpath = os.path.expanduser('~')
        
        self.filepath = self.userpath #os.getcwd()
        self.systemfile = None
        self.plotfile = None
        self.resultfile = None
        
        self.numprocessors = max(2,multiprocessing.cpu_count())
        self.path_base = os.path.dirname(os.path.realpath(__file__))
        self.path_workingdir = os.path.join(self.path_base,"../calc/")
        self.path_cpp_real = os.path.join(self.path_base,"../calc/pairinteraction-real")
        self.path_cpp_complex = os.path.join(self.path_base,"../calc/pairinteraction-complex")
        
        if os.name == 'nt': self.path_out = os.path.join(self.userpath, "pairinteraction/")
        else: self.path_out = os.path.join(self.userpath, ".pairinteraction/")
        self.path_cache = os.path.join(self.path_out, "cache/")
        self.path_lastsettings = os.path.join(self.path_out, "lastsettings/")
        
        self.path_system_last = os.path.join(self.path_lastsettings,"lastsettings.sconf")
        self.path_plot_last = os.path.join(self.path_lastsettings,"lastsettings.pconf")
        self.path_view_last = os.path.join(self.path_lastsettings,"lastsettings.json")
        
        self.path_config = os.path.join(self.path_out,"conf.json")
        self.path_version = os.path.join(self.path_out,"version.json")
        
        self.proc = None
        
        self.thread = Worker()
        self.timer = QtCore.QTimer()
        
        self.momentumcolors = [(55,126,184),(77,175,74),(228,26,28),(152,78,163),(0,0,0),(255//5,255//5,255//5)] # s, p, d, f, other, undetermined
        self.symmetrycolors = [(0,0,0),(140,81,10),(1,102,94)] # all, sym, asym
        
        self.momentummat = [None]*3
        self.labelmat = [None]*3
        self.labelstates = [None]*3
        self.momentumstrings = [None]*3
        self.stateidx_field = [None]*3
        self.yMin_field = [None]*3
        self.yMax_field = [None]*3
        
        self.buffer_basis = [{},{},{}]
        self.buffer_energies = [{},{},{}]
        self.buffer_positions = [{},{},{}]
        self.buffer_boolarr = [{},{},{}]
        self.buffer_basis_potential = {}
        self.buffer_energies_potential = {}
        self.buffer_positions_potential = {}
        
        self.buffer_energiesMap = [{},{},{}]
        self.buffer_positionsMap = [{},{},{}]
        self.buffer_overlapMap = [{},{},{}]
        self.buffer_energiesMap_potential = {}
        self.buffer_positionsMap_potential = {}
        self.buffer_overlapMap_potential = {}
        
        self.lines_buffer_minIdx = {}
        self.colormap_buffer_minIdx_potential = 0
        self.colormap_buffer_minIdx_field = [0]*3
        self.lines_buffer_minIdx_field = [0]*3
        

        
        
        #clrmp = pg.ColorMap(pos,color)
        #self.lut = clrmp.getLookupTable()
                
        self.ui.gradientwidget_plot_gradient.setOrientation("top")
        self.ui.gradientwidget_plot_gradient.loadPreset('cubehelix1')
        
        self.tab_field2 = self.ui.tabwidget_plotter.widget(1)
        
        
        self.storage_data = [[],[],[]]
        self.storage_states = [None, None, None]
        self.storage_configuration = [[None, None], [None, None], [None, None]]

        # TODOs
        self.ui.lineedit_system_theta.setEnabled(False)
        self.ui.lineedit_system_precision.setEnabled(False)
        self.ui.checkbox_system_dq.setEnabled(False)
        self.ui.checkbox_system_qq.setEnabled(False)
        self.ui.checkbox_system_qq.setEnabled(False)
        self.ui.combobox_plot_overlapTransition.setEnabled(False)
        self.ui.radiobutton_plot_overlapTransition.setEnabled(False)
        
        
        
                
        

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
        self.ui.checkbox_system_samebasis.toggled.connect(self.toggleSamebasis)
        
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
        
        self.ui.pushbutton_field1_save.clicked.connect(self.saveResult)
        self.ui.pushbutton_field2_save.clicked.connect(self.saveResult)
        self.ui.pushbutton_potential_save.clicked.connect(self.saveResult)
        
        self.timer.timeout.connect(self.checkForData)
        
        
        # Check version
        
        # Load version
        version_settings_saved = None
        version_cache_saved = None
        if os.path.isfile(self.path_version):
            with open(self.path_version, 'r') as f:
                params = json.load(f)
                version_settings_saved = params["version_settings"]
                version_cache_saved =  params["version_cache"]
        
        # Compare version
        if os.path.exists(self.path_out) and version_settings_saved != version_settings and version_cache_saved != version_cache:
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
                sys.exit()
        
        elif os.path.exists(self.path_lastsettings) and version_settings_saved != version_settings:
            msg = QtGui.QMessageBox()
            msg.setText('A new program version has been installed. Due to configuration changes, settings have to be cleared. This deletes the directory {}.'.format(self.path_lastsettings))
            msg.setIcon(QtGui.QMessageBox.Information);
            msg.addButton(QtGui.QMessageBox.Cancel)
            msg.addButton(QtGui.QMessageBox.Ok)
            msg.setDefaultButton(QtGui.QMessageBox.Ok)
            answer = msg.exec()
                
            # Delete directory
            if answer == QtGui.QMessageBox.Ok:
                shutil.rmtree(self.path_lastsettings)
            else:
                sys.exit()
        
        elif os.path.exists(self.path_cache) and version_cache_saved != version_cache:
            msg = QtGui.QMessageBox()
            msg.setText('A new program version has been installed. Due to database changes, the cache has to be cleared. This deletes the directory {}.'.format(self.path_cache))
            msg.setIcon(QtGui.QMessageBox.Information);
            msg.addButton(QtGui.QMessageBox.Cancel)
            msg.addButton(QtGui.QMessageBox.Ok)
            msg.setDefaultButton(QtGui.QMessageBox.Ok)
            answer = msg.exec()
                
            # Delete directory
            if answer == QtGui.QMessageBox.Ok:
                shutil.rmtree(self.path_cache)
            else:
                sys.exit()
        
        '''msg = None
        
        if os.path.exists(self.path_lastsettings):
            # Load version
            version_settings_saved = None
            if os.path.isfile(self.path_version):
                with open(self.path_version, 'r') as f:
                    version_settings_saved = json.load(f)["version_settings"]
            
            # Compare version
            if version_settings_saved != version_settings:
                """msg = QtGui.QMessageBox() # TODO make class
                msg.setText('A new version of the program has been installed. Clear cache and user data to avoid compatibility issues? This deletes the directory "{}".'.format(self.path_out))
                msg.setIcon(QtGui.QMessageBox.Information);
                checkbox = QtGui.QCheckBox("don't show message again")
                checkbox.blockSignals(True)
                msg.addButton(checkbox, QtGui.QMessageBox.ApplyRole)
                msg.addButton(QtGui.QMessageBox.Yes)
                msg.addButton(QtGui.QMessageBox.No)
                msg.setDefaultButton(QtGui.QMessageBox.Yes)
                answer = msg.exec()
                
                # Save version if "don't show message again")
                if answer == QtGui.QMessageBox.No and checkbox.isChecked():
                    with open(self.path_version, 'w') as f:
                        json.dump({'version_settings': version_settings, 'version_cache': version_cache}, f, indent=4, sort_keys=True)
                
                # Delete directory
                if answer == QtGui.QMessageBox.Yes:
                    shutil.rmtree(self.path_out)"""
                
                msg = QtGui.QMessageBox() # TODO make class
                #msg.setText('A new program version has been installed. Due to major changes, cache and last settings have to be cleared. This deletes the directory "{}".'.format(self.path_out))
                msg.setText('A new program version has been installed. Due to configuration changes, the last settings have to be cleared. This deletes the directory "{}".'.format(self.path_lastsettings))
                msg.setIcon(QtGui.QMessageBox.Information);
                msg.addButton(QtGui.QMessageBox.Cancel)
                msg.addButton(QtGui.QMessageBox.Ok)
                msg.setDefaultButton(QtGui.QMessageBox.Ok)
                answer = msg.exec()
                    
                # Delete directory
                if answer == QtGui.QMessageBox.Ok:
                    shutil.rmtree(self.path_lastsettings)
                else:
                    sys.exit()
                    
        if os.path.exists(self.path_cache):     
            # Load version
            version_cache_saved = None
            if os.path.isfile(self.path_version):
                with open(self.path_version, 'r') as f:
                    version_cache_saved = json.load(f)["version_cache"]
                        
            # Compare version
            if version_cache_saved != version_cache:
                msg = QtGui.QMessageBox()
                msg.setText('A new program version has been installed. Due to database changes, the cache has to be cleared. This deletes the directory "{}".'.format(self.path_cache))
                msg.setIcon(QtGui.QMessageBox.Information);
                msg.addButton(QtGui.QMessageBox.Cancel)
                msg.addButton(QtGui.QMessageBox.Ok)
                msg.setDefaultButton(QtGui.QMessageBox.Ok)
                answer = msg.exec()
                
                # Delete directory
                if answer == QtGui.QMessageBox.Ok:
                    shutil.rmtree(self.path_cache)
                else:
                    sys.exit()'''
            
        # Create directories
        if not os.path.exists(self.path_out):
            os.makedirs(self.path_out)	
            if os.name == 'nt':
                ret = ctypes.windll.kernel32.SetFileAttributesW(self.path_out,FILE_ATTRIBUTE_HIDDEN)
                if not ret: raise ctypes.WinError()
            
            with open(self.path_version, 'w') as f:
                json.dump({'version_settings': version_settings, 'version_cache': version_cache}, f, indent=4, sort_keys=True)
        
        if not os.path.exists(self.path_lastsettings):
            os.makedirs(self.path_lastsettings)
        
            with open(self.path_version, 'r') as f:
                version_cache_saved = json.load(f)["version_cache"]
            
            with open(self.path_version, 'w') as f:
                json.dump({'version_settings': version_settings, 'version_cache': version_cache_saved}, f, indent=4, sort_keys=True)
            
        if not os.path.exists(self.path_cache):
            os.makedirs(self.path_cache)
        
            with open(self.path_version, 'r') as f:
                version_settings_saved = json.load(f)["version_settings"]
            
            with open(self.path_version, 'w') as f:
                json.dump({'version_settings': version_settings_saved, 'version_cache': version_cache}, f, indent=4, sort_keys=True)
        
        # Load last settings
        if not os.path.isfile(self.path_system_last):
            copyfile(os.path.join(self.path_base,"example.sconf"),self.path_system_last)
        with open(self.path_system_last, 'r') as f:
            params = json.load(f)
            self.systemdict.load(params)
    
        if not os.path.isfile(self.path_plot_last):
            copyfile(os.path.join(self.path_base,"example.pconf"),self.path_plot_last)
        with open(self.path_plot_last, 'r') as f:
            params = json.load(f)
            self.ui.gradientwidget_plot_gradient.restoreState(params["gradientwidget"])
            del params["gradientwidget"]
            self.plotdict.load(params)
    
        if os.path.isfile(self.path_view_last):
            with open(self.path_view_last, 'r') as f:
                params = json.load(f)
                self.ui.tabwidget_config.setCurrentIndex(params["config"])
                self.ui.tabwidget_plotter.setCurrentIndex(params["plotter"])
                self.ui.toolbox_system.setCurrentIndex(params["system"])
                if "filepath" in params.keys(): self.filepath = params["filepath"]
        
        
        
        
        self.samebasis_state = self.ui.checkbox_system_samebasis.checkState()
                
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
        self.ui.checkbox_system_samebasis.toggled.emit(self.ui.checkbox_system_samebasis.isChecked())
        
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
            plotarea.scene().contextMenu = None
            plotarea.plotItem.ctrlMenu = None
            #plotarea.getViewBox().menu = None # AttributeError: 'NoneType' object has no attribute 'popup'
            
            plotarea.getAxis("bottom").setZValue(1000) # HACK to bring axis into the foreground
            plotarea.getAxis("left").setZValue(1000) # HACK to bring axis into the foreground
        
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
            
    """def resizeEvent(self, event):
        super().resizeEvent(event)
        if sys.platform == "darwin": QtGui.QApplication.processEvents() # hack to circumvent the no-redraw-after-resizing-bug"""
    
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
        
        # === process field and potential maps ===
             
        for idx in range(3):
            
            if idx > 0 and not self.thread.dataqueue_field1.empty(): continue
            if idx > 1 and not self.thread.dataqueue_field2.empty(): continue
            
            basisfile = [self.thread.basisfile_field1, self.thread.basisfile_field2, self.thread.basisfile_potential][idx]
            dataqueue = [self.thread.dataqueue_field1, self.thread.dataqueue_field2, self.thread.dataqueue_potential][idx]
               
            # --- load basis states ---
                        
            if basisfile != "":
                # load basis
                basis = np.loadtxt(basisfile)
                
                # save basis
                self.storage_states[idx] = basis
                
                # determine which state to highlite
                if self.ui.groupbox_plot_overlap.isChecked():
                    if idx == 0:
                        try: self.stateidx_field[0] = np.where(np.all(basis[:,[1,2,3,4]] == self.overlapstate[None,[0,1,2,3]],axis=-1))[0][0]
                        except: self.stateidx_field[0] = -1
                    if  idx == 1 or (idx == 0 and self.thread.samebasis):
                        try: self.stateidx_field[1] = np.where(np.all(basis[:,[1,2,3,4]] == self.overlapstate[None,[4,5,6,7]],axis=-1))[0][0]
                        except: self.stateidx_field[1] = -1
                    elif idx == 2:
                        try: self.stateidx_field[2] = np.where(np.all(basis[:,[1,2,3,4,5,6,7,8]] == self.overlapstate[None,:],axis=-1))[0][0]
                        except: self.stateidx_field[2] = -1
                
                    self.yMin_field[idx] = None
                    self.yMax_field[idx] = None
            
                # calculate a matrix that can be used to determine the momenta inside a basis element 
                if idx == 0 or idx == 1:
                    momentum = basis[:,2]
                    self.momentummat[idx] = sparse.csc_matrix((momentum[:,None] == np.arange(np.max(momentum)+1)[None,:]).astype(int))

                # extract labels
                # TODO only i necessary !!!!!
                
                if idx == 0 or idx == 1:
                    nlj = basis[:,[1,2,3]]
                elif idx == 2:
                    nlj = basis[:,[1,2,3,5,6,7]]
                
                # sort pair state names
                if idx == 2 and self.thread.samebasis:
                    firstsmaller = np.argmax(np.append(nlj[:,0:3] < nlj[:,3:6],np.ones((len(nlj),1),dtype=bool),axis=-1),axis=-1) # TODO in Funktion auslagern
                    firstsmaller_reverted = np.argmax(np.append(nlj[:,3:6] < nlj[:,0:3],np.ones((len(nlj),1),dtype=bool),axis=-1),axis=-1) # TODO in Funktion auslagern
                    namesToSwap = firstsmaller > firstsmaller_reverted
                    nlj[namesToSwap] = nlj[namesToSwap][:,[3,4,5,0,1,2]]
                
                sorter = np.lexsort(nlj.T[::-1])
                nlj = nlj[sorter]
                diff = np.append([True],np.diff(nlj, axis=0).any(axis=1))
                cumsum = np.cumsum(diff)[np.argsort(sorter)]
            
                # determine labels of states
                self.labelstates[idx] = nlj[diff]            
                self.labelmat[idx] = sparse.coo_matrix((np.ones_like(cumsum),(np.arange(len(cumsum)),cumsum-1)),shape=(len(cumsum), len(self.labelstates[idx]))).tocsr() #nStates, nLabels
                
                """try:
                    if idx == 0: self.labelidx1 = np.where(np.all(basis[:,[1,2,3,4]] == self.unperturbedstate[None,[0,1,2,3]],axis=-1))[0][0]
                    elif idx == 1: self.labelidx1 = np.where(np.all(basis[:,[1,2,3,4]] == self.unperturbedstate[None,[4,5,6,7]],axis=-1))[0][0]
                    elif idx == 2: self.labelidx1 = np.where(np.all(basis[:,[1,2,3,4,5,6,7,8]] == self.unperturbedstate[None,:],axis=-1))[0][0]
                except:
                    self.labelidx1 =  -1
                
                try:
                    if idx == 0 and self.thread.samebasis: self.labelidx2 = np.where(np.all(basis[:,[1,2,3,4]] == self.unperturbedstate[None,[4,5,6,7]],axis=-1))[0][0]
                    else: raise
                except:
                    self.labelidx2 = -1"""
 
                # determine labels of momenta
                if idx == 0 or idx == 1: # TODO !!!
                    self.momentumstrings[idx] = [" {}".format(i) for i in np.arange(np.max(self.labelstates[idx][:,1])+1).astype(np.int)]
                elif idx == 2:
                    self.momentumstrings[idx] = [" {}".format(i) for i in np.arange(np.max(self.labelstates[idx][:,[1,4]])+1).astype(np.int)]
                self.momentumstrings[idx][:4] = ['S','P','D','F']
            
                # remove basis file from hard disk
                os.remove(basisfile)
                
                # clear variables
                self.lines_buffer_minIdx_field = {}
                self.buffer_basis = {}
                self.buffer_energies = {}
                self.buffer_positions = {}
                self.buffer_boolarr = {}
            
                # indicate that the basis file is already processed
                if idx == 0: self.thread.basisfile_field1 = ""
                elif idx == 1: self.thread.basisfile_field2 = ""
                elif idx == 2: self.thread.basisfile_potential = ""
            
            # --- check if there is some new data and if yes, plot it ---
            
            if not dataqueue.empty():
                
                graphicsview_plot = [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot, self.ui.graphicsview_potential_plot]
                
                numBlocks = [self.thread.numBlocks_field1, self.thread.numBlocks_field2, self.thread.numBlocks_potential][idx]
                minE = [self.minE_field1, self.minE_field2, self.minE_potential][idx]
                maxE = [self.maxE_field1, self.maxE_field2, self.maxE_potential][idx]
                
                # --- storage that allows to draw the hole buffer at once, at least if it is no very large ---
                x = np.array([])
                y = np.array([])
                l = np.array([])
                s = np.array([])
            
                while not dataqueue.empty() and dataamount < 5000: # stop loop if enough data is collected
                
                    # --- load eigenvalues (energies, y value) and eigenvectors (basis) ---
                    filestep, blocknumber, filename = dataqueue.get()
                    
                    # save data
                    self.storage_data[idx].append([filestep, blocknumber, filename])
                    
                    eigensystem = Eigensystem(filename)
                    energies = eigensystem.energies
                    basis = eigensystem.basis
                    
                    if idx == 2: symmetry = blocknumber % 3
                                    
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
                    if idx == 0 or idx == 1:
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
                    if idx == 0 or idx == 1:
                        if self.constEField and not self.constBField:
                            position = self.get1DPosition(float(eigensystem.params["Bx"]),float(eigensystem.params["By"]),float(eigensystem.params["Bz"]))*self.converter_x_field
                        else:
                            position = self.get1DPosition(float(eigensystem.params["Ex"]),float(eigensystem.params["Ey"]),float(eigensystem.params["Ez"]))*self.converter_x_field
                    elif idx == 2:
                        if self.constDistance and not self.constEField:
                            position = self.get1DPosition(float(eigensystem.params["Ex"]),float(eigensystem.params["Ey"]),float(eigensystem.params["Ez"]))*self.converter_x_potential
                        elif self.constDistance and not self.constBField:
                            position = self.get1DPosition(float(eigensystem.params["Bx"]),float(eigensystem.params["By"]),float(eigensystem.params["Bz"]))*self.converter_x_potential
                        else:
                            position = float(eigensystem.params["R"])*self.converter_x_potential

                    # --- draw labels at the beginning of the plotting ---
                    if self.ui.groupbox_plot_labels.isChecked() and filestep == 0:
                
                        # probability to find a label inside a basis element
                        if not hasattr(self, 'labelprob_energy') or self.labelprob_energy is None:
                            self.labelprob = (probs.T*self.labelmat[idx]).tocsr() #nBasis, nLabels # TODO !!! tocsr
                            self.labelprob_energy = [energies]
                        else:
                            csr_vappend(self.labelprob,(probs.T*self.labelmat[idx]).tocsr())
                            self.labelprob_energy.append(energies)
                        
                        labelprob_num_potential = len(self.labelprob_energy)
                                                                        
                        if labelprob_num_potential >= numBlocks:
                            
                            self.labelprob_energy = np.concatenate(self.labelprob_energy)
                        
                            # total probability to find a label
                            cumprob = np.array(self.labelprob.sum(axis=0).flat)
                            boolarr = cumprob > 0.1
                                                         
                            # normalize in such a way that the total probability is always one
                            idxarr, = np.nonzero(boolarr)
                            normalizer = sparse.coo_matrix((1/cumprob[idxarr],(idxarr,np.arange(len(idxarr)))),shape=(self.labelprob.shape[1], len(idxarr))).tocsr()
                        
                            # store the energetic expectation value of the labels
                            labelenergies = ((self.labelprob*normalizer).T*self.labelprob_energy)
                            
                            if len(labelenergies) == 0: continue
                                                    
                            # store the position of the labels
                            labelposition = position
    
                            # get size and alpha value of labels
                            size = '{}'.format(max(int(round(self.ui.spinbox_plot_szLabel.value()*11)),1))
                            alpha = int(round(self.ui.spinbox_plot_transpLabel.value()*255))
                        
                            # draw the labels
                            for labelstate, labelenergy in zip(self.labelstates[idx][boolarr],labelenergies):
                                
                                if ((idx == 0 or idx == 1) and self.leftSmallerRight) or (idx == 2 and self.leftSmallerRight_potential):
                                    anchorX = 0
                                else:
                                    anchorX = 1
                                
                                if (idx == 0 and np.all(labelstate == self.unperturbedstate[[0,1,2]])) \
                                        or ((idx == 1 or self.thread.samebasis) and np.all(labelstate == self.unperturbedstate[[4,5,6]])) \
                                        or (idx == 2 and np.all(labelstate == self.unperturbedstate[[0,1,2,4,5,6]])) \
                                        or ((idx == 2 and self.thread.samebasis) and np.all(labelstate == self.unperturbedstate[[4,5,6,0,1,2]])):
                                    color_fill = (255,192,203,alpha)
                                    color_border = (255,182,193,255)
                                    zvalue = 16
                                else:
                                    color_fill = (250,235,215,alpha)
                                    color_border = (255,228,181,255)
                                    zvalue = 15
                                    
                                if idx == 0 or idx == 1:
                                    sn, sl, sj = labelstate
                                    text = pg.TextItem(html='<div style="text-align: center; font-size: '+size+'pt;">'\
                                        +'<span style="color: rgba(0,0,0,255);">{}{}<sub style="font-size: '.format(int(sn),self.momentumstrings[idx][int(sl)]) \
                                        +size+'pt;">{}/2</sub></span></div>'.format(int(2*sj)),anchor=(anchorX, 0.5),fill=color_fill,border=color_border)
                                elif idx == 2:
                                    sn1, sl1, sj1, sn2, sl2, sj2 = labelstate
                                    text = pg.TextItem(html='<div style="text-align: center; font-size: '+size+'pt;"><span style="color: rgba(0,0,0,255);">'\
                                        +'{}{}<sub style="font-size: '.format(int(sn1),self.momentumstrings[idx][int(sl1)]) \
                                        +size+'pt;">{}/2</sub>'.format(int(2*sj1))\
                                        +' {}{}<sub style="font-size: '.format(int(sn2),self.momentumstrings[idx][int(sl2)]) \
                                        +size+'pt;">{}/2</sub>'.format(int(2*sj2))\
                                        +'</span></div>',anchor=(anchorX, 0.5),fill=color_fill,border=color_border)
                                
                                text.setPos(labelposition, labelenergy)
                                text.setZValue(zvalue)
                                graphicsview_plot[idx].addItem(text)
                                
                            posx = labelposition*np.ones_like(labelenergies)+1e-12
                            posy = labelenergies+1e-12
                            curve = PointsItem(np.append(posx, posx-2e-12), np.append(posy, posy-2e-12), 0, 0, (255,255,255))
                            curve.setZValue(5)
                            graphicsview_plot[idx].addItem(curve)               
                        
                            # drawing labels take some time
                            dataamount += 3000
                            
                            self.labelprob = None
                            self.labelprob_energy = None
                    
                    # --- draw color map ---
                    if self.ui.groupbox_plot_overlap.isChecked():
                        # --- get settings ---
                        # get size and alpha value
                        size = self.ui.spinbox_plot_szOverlap.value()
                        alpha = min(self.ui.spinbox_plot_transpOverlap.value(),0.9999) # HACK
                        
                        # get resolution
                        res = self.ui.spinbox_plot_resolution.value()
                        
                        # calculate size of color map
                        height_pixelunits = res
                        enlargement = int(max(np.round((height_pixelunits/self.steps-2)/2),0))
                        width_pixelunits = 5+4*enlargement
                        
                        # calculate values to smooth the colormap
                        smootherX = (enlargement*2+2)*1/2
                        smootherY = height_pixelunits*1/150*size
                        
                        # --- build buffers ---
                        # initialize arrays if first run at a new position ("filestep")
                        if filestep not in self.buffer_positionsMap[idx].keys():
                            self.buffer_positionsMap[idx][filestep] = position
                            self.buffer_energiesMap[idx][filestep] = []
                            self.buffer_overlapMap[idx][filestep] = []
                            
                        # try to get limits
                        if self.yMax_field[idx] is None and maxE is not None:
                            self.yMax_field[idx]  = maxE*self.converter_y  
                        if self.yMin_field[idx] is None and minE is not None:
                            self.yMin_field[idx] = minE*self.converter_y
                                                                    
                        # check if limits do not exists
                        if self.yMax_field[idx] is None or self.yMin_field[idx] is None:
                            # append the energies to the arrays
                            self.buffer_energiesMap[idx][filestep].append(energies)
                            
                            # append the overlaps to the arrays 
                            if self.stateidx_field[idx] >= 0:
                                self.buffer_overlapMap[idx][filestep].append(probs[self.stateidx_field[idx]].toarray().flatten())
                            else:
                                self.buffer_overlapMap[idx][filestep].append(np.zeros_like(energies))
                                
                            if self.thread.samebasis and self.stateidx_field[1] >= 0 and self.stateidx_field[1] != self.stateidx_field[0]:
                                self.buffer_overlapMap[idx][filestep][-1] += probs[self.stateidx_field[1]].toarray().flatten()
                                
                            # check if all data of the zeroth position is collected
                            if 0 in self.buffer_overlapMap[idx].keys() and len(self.buffer_overlapMap[idx][0]) == numBlocks:
                                # make limits
                                if self.yMin_field[idx] is None:
                                    self.yMin_field[idx] = np.nanmin(np.concatenate(self.buffer_energiesMap[idx][0]))
                                if self.yMax_field[idx] is None:
                                    self.yMax_field[idx] = np.nanmax(np.concatenate(self.buffer_energiesMap[idx][0]))
                                
                                # calculate energy-indices
                                for f in self.buffer_energiesMap[idx].keys():
                                    for i in range(len(self.buffer_energiesMap[idx][f])):
                                        boolarr = (self.buffer_energiesMap[idx][f][i] >= self.yMin_field[idx]) & (self.buffer_energiesMap[idx][f][i] <= self.yMax_field[idx])
                                        self.buffer_energiesMap[idx][f][i] = bytescale(self.buffer_energiesMap[idx][f][i][boolarr], low=0, high=height_pixelunits-1, cmin=self.yMin_field[idx], cmax=self.yMax_field[idx])
                                        self.buffer_overlapMap[idx][f][i] = self.buffer_overlapMap[idx][f][i][boolarr]
                        else:
                            # determine whether the energies lie within the limits
                            boolarr = (energies >= self.yMin_field[idx]) & (energies <= self.yMax_field[idx])
    
                            # append the energy-indices to the arrays
                            self.buffer_energiesMap[idx][filestep].append(bytescale(energies[boolarr], low=0, high=height_pixelunits-1, cmin=self.yMin_field[idx], cmax=self.yMax_field[idx]))
                            
                            # append the overlaps to the arrays
                            if self.stateidx_field[idx] >= 0:
                                self.buffer_overlapMap[idx][filestep].append(probs[self.stateidx_field[idx]].toarray().flatten()[boolarr])
                            else:
                                self.buffer_overlapMap[idx][filestep].append(np.zeros_like(energies[boolarr]))
                            
                            if self.thread.samebasis and self.stateidx_field[1] >= 0 and self.stateidx_field[1] != self.stateidx_field[0]:
                                self.buffer_overlapMap[idx][filestep][-1] += probs[self.stateidx_field[1]].toarray().flatten()[boolarr]
                        
                        # --- build color maps starting at the lowest position---                    
                        # loop over positions ("filestep") as long as three subsequent positions ("self.colormap_buffer_minIdx_field[idx]+0,1,2") are within the buffers
                        while  True:
                            # break if limits do not exist
                            if self.yMax_field[idx] is None or self.yMin_field[idx] is None: break
                            
                            # break if buffer index is not in the buffer
                            bufferidx = []
                            
                            if self.colormap_buffer_minIdx_field[idx] != 0: # not start
                                if self.colormap_buffer_minIdx_field[idx]-1 not in self.buffer_positionsMap[idx].keys(): break
                                bufferidx.append(-1)
                                
                            if self.colormap_buffer_minIdx_field[idx] not in self.buffer_positionsMap[idx].keys(): break
                            bufferidx.append(0)
                                                    
                            if self.colormap_buffer_minIdx_field[idx] != self.steps-1: # not end
                                if self.colormap_buffer_minIdx_field[idx]+1 not in self.buffer_positionsMap[idx].keys(): break
                                bufferidx.append(1)
                            
                            # break if the data is not buffered of all blocks, yet
                            tobreak = False
                            for i in bufferidx:
                                if len(self.buffer_energiesMap[idx][self.colormap_buffer_minIdx_field[idx]+i]) < numBlocks: tobreak = True
                            if tobreak: break
 
                            # calculate position-indices
                            positions = [self.buffer_positionsMap[idx][self.colormap_buffer_minIdx_field[idx]+i] for i in bufferidx]
                            
                            # add outer points if at the end or start
                            if self.colormap_buffer_minIdx_field[idx] == self.steps-1: # end
                                positions = positions + [2*positions[1]-positions[0]]
                            elif self.colormap_buffer_minIdx_field[idx] == 0: # start
                                positions = [2*positions[0]-positions[1]] + positions
                                
                            # determine limits of the color map part
                            posLeft = (positions[0]+positions[1])/2
                            posRight = (positions[1]+positions[2])/2
                            idx_left, idx_right = bytescale(np.array([posLeft,posRight]), low=0, high=width_pixelunits-1, cmin=positions[0], cmax=positions[-1])
                            
                            # calculate unit converters
                            self.displayunits2pixelunits_x = width_pixelunits/(positions[2]-positions[0])
                            self.displayunits2pixelunits_y = height_pixelunits/(self.yMax_field[idx]-self.yMin_field[idx])
                            
                            # build map
                            colormap = np.zeros((width_pixelunits,height_pixelunits)) # x-y-coordinate system, origin is at the bottom left corner 
                            
                            for i in bufferidx:
                                overlap = np.concatenate(self.buffer_overlapMap[idx][self.colormap_buffer_minIdx_field[idx]+i])
                                pos = self.buffer_positionsMap[idx][self.colormap_buffer_minIdx_field[idx]+i]
                                idx_pos = bytescale(pos, low=0, high=width_pixelunits-1, cmin=positions[0], cmax=positions[-1])
                                idx_energies = np.concatenate(self.buffer_energiesMap[idx][self.colormap_buffer_minIdx_field[idx]+i])
                                
                                if self.logscale:
                                    overlap[overlap < 1e-2] = 1e-2
                                    overlap = (2+np.log10(overlap))/2
                                
                                colormap[idx_pos,:] = sparse.coo_matrix((overlap,(idx_energies,np.arange(len(idx_energies)))),shape=(height_pixelunits, len(idx_energies))).sum(axis=-1).flat
                                
                                dataamount += len(idx_energies)*10
    
                            # smoothing
                            colormap = gaussian_filter(colormap,(smootherX,smootherY),mode='constant')
                            
                            # cutting
                            colormap = colormap[idx_left:idx_right]
                            
                            # normalizing
                            normalizer = np.zeros((width_pixelunits,int(2*3*smootherY+1)))
                            normalizer[int((width_pixelunits-1)/2),int(3*smootherY)] = 1
                            normalizer = gaussian_filter(normalizer,(smootherX,smootherY),mode='constant')
                            
                            colormap /= np.max(normalizer)
                            
                            # plotting
                            img = pg.ImageItem(image=colormap, opacity=alpha, autoDownsample=True, lut=self.lut, levels=[-0.002,1]) # HACK
                            img.setRect(QtCore.QRectF(posLeft-0.5/self.displayunits2pixelunits_x,self.yMin_field[idx]-0.5/self.displayunits2pixelunits_y,posRight-posLeft,self.yMax_field[idx]-self.yMin_field[idx]+1/self.displayunits2pixelunits_y)) # xMin, yMin_field[idx], xSize, ySize # TODO energyMin anpassen wegen Pixelgroesse
                            img.setZValue(3)
                            graphicsview_plot[idx].addItem(img)
                                                            
                            # remove plotted data from buffer
                            if self.colormap_buffer_minIdx_field[idx] != 0: # not start
                                del self.buffer_energiesMap[idx][self.colormap_buffer_minIdx_field[idx]-1]
                                del self.buffer_positionsMap[idx][self.colormap_buffer_minIdx_field[idx]-1]
                                del self.buffer_overlapMap[idx][self.colormap_buffer_minIdx_field[idx]-1]
                            
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
                    
                        self.buffer_basis[blocknumber][filestep] = basis
                        self.buffer_energies[blocknumber][filestep] = energies
                        self.buffer_positions[blocknumber][filestep] = position
                        self.buffer_boolarr[blocknumber][filestep] = []
                        
                        if idx == 0 or idx == 1:
                            # cut the momenta to reasonable values
                            momentum[momentum>len(self.momentumcolors)-2] = len(self.momentumcolors)-2
                            momentum[momentum<0] = len(self.momentumcolors)-1
                        
                            # loop over momenta
                            for i in range(len(self.momentumcolors)):
                                # determine which basis elements have the current momentum
                                boolarr = momentum == i
                                self.buffer_boolarr[blocknumber][filestep].append(boolarr)
                        elif idx == 2:
                            self.buffer_boolarr[blocknumber][filestep].append(np.ones_like(self.buffer_energies[blocknumber][filestep], dtype = np.bool))
                                             
                        # get size and alpha value of points
                        size = self.ui.spinbox_plot_szLine.value()
                        alpha = self.ui.spinbox_plot_transpLine.value()*255
                        
                        """# legend
                        if idx == 2 and filestep == 0 and symmetry != 0:
                            graphicsview_plot[idx].addLegend()
                            style = pg.PlotDataItem(pen = pg.mkPen(self.symmetrycolors[1]+(alpha,),width=size,cosmetic=True))
                            graphicsview_plot[idx].plotItem.legend.addItem(style, "sym")
                            style = pg.PlotDataItem(pen = pg.mkPen(self.symmetrycolors[2]+(alpha,),width=size,cosmetic=True))
                            graphicsview_plot[idx].plotItem.legend.addItem(style, "asym")"""
                                                                                                
                        while self.lines_buffer_minIdx_field[blocknumber] in self.buffer_basis[blocknumber].keys() and self.lines_buffer_minIdx_field[blocknumber]+1 in self.buffer_basis[blocknumber].keys():
                            # determine the data to plot
                            overlap = np.abs(self.buffer_basis[blocknumber][self.lines_buffer_minIdx_field[blocknumber]].conj().T*self.buffer_basis[blocknumber][self.lines_buffer_minIdx_field[blocknumber]+1]) # nBasis first, nBasis second
                            
                            overlap.data[overlap.data <= np.sqrt(self.ui.spinbox_plot_connectionthreshold.value())] = 0
                            overlap.eliminate_zeros()
                            csr_keepmax(overlap)
                            
                            
                            overlap = overlap.tocoo()

                            iFirst = overlap.row
                            iSecond = overlap.col
                            
                            ydata = np.transpose([self.buffer_energies[blocknumber][self.lines_buffer_minIdx_field[blocknumber]][iFirst],self.buffer_energies[blocknumber][self.lines_buffer_minIdx_field[blocknumber]+1][iSecond]])
                            xdata = np.ones_like(ydata)
                            xdata[:,0] *= self.buffer_positions[blocknumber][self.lines_buffer_minIdx_field[blocknumber]]
                            xdata[:,1] *= self.buffer_positions[blocknumber][self.lines_buffer_minIdx_field[blocknumber]+1]
                        
                            # loop over momenta
                            numMomenta = len(self.buffer_boolarr[blocknumber][self.lines_buffer_minIdx_field[blocknumber]])
                            
                            for i in range(numMomenta):
                                boolarr = self.buffer_boolarr[blocknumber][self.lines_buffer_minIdx_field[blocknumber]][i][iFirst]
                                if np.sum(boolarr) == 0: continue
                            
                                # determine the associated color
                                if idx == 0 or idx == 1:
                                    color = self.momentumcolors[i]
                                elif idx == 2:
                                    color = self.symmetrycolors[symmetry]
                                
                                # plot the data
                                curve = MultiLine(xdata[boolarr], ydata[boolarr], size, alpha, color) # TODO alpha and color der Funktion zusammen uebergeben
                                curve.setZValue(7)
                                graphicsview_plot[idx].addItem(curve)
                            
                            del self.buffer_basis[blocknumber][self.lines_buffer_minIdx_field[blocknumber]]
                            del self.buffer_energies[blocknumber][self.lines_buffer_minIdx_field[blocknumber]]
                            del self.buffer_positions[blocknumber][self.lines_buffer_minIdx_field[blocknumber]]
                            del self.buffer_boolarr[blocknumber][self.lines_buffer_minIdx_field[blocknumber]]
                            
                            # increase the buffer index
                            self.lines_buffer_minIdx_field[blocknumber] += 1
                            
                            dataamount += len(iFirst)*10
                    
                    # --- store data to plot several points at once ---
                    if self.ui.groupbox_plot_points.isChecked():
                        x = np.append(x,position*np.ones_like(energies))
                        y = np.append(y,energies)
                        if idx == 0 or idx == 1:
                            l = np.append(l,momentum)
                        elif idx == 2:
                            s = np.append(s,symmetry*np.ones_like(energies))
                        
                        dataamount += len(x)
                    
                # --- plot the stored data ---
                if self.ui.groupbox_plot_points.isChecked() and len(x) > 0:
            
                    # get size and alpha value of points
                    size = self.ui.spinbox_plot_szPoint.value()
                    alpha = self.ui.spinbox_plot_transpPoint.value()*255
                    
                    if idx == 0 or idx == 1:
                        # cut the momenta to reasonable values
                        l[l>len(self.momentumcolors)-2] = len(self.momentumcolors)-2
                        l[l<0] = len(self.momentumcolors)-1
                
                    if idx == 0 or idx == 1:
                        looprange = len(self.momentumcolors)
                    elif idx == 2:
                        looprange = len(self.symmetrycolors)
                
                    # loop over momenta
                    for i in range(looprange):
                        if idx == 0 or idx == 1:
                            # determine which basis elements have the current momentum
                            boolarr = l == i
                            if (np.sum(boolarr) == 0): continue
                    
                            # determine the associated color
                            color = self.momentumcolors[i]
                        elif idx == 2:
                            # determine which basis elements have the current symmetry
                            boolarr = s == i
                            if (np.sum(boolarr) == 0): continue
                    
                            # determine the associated color
                            color = self.symmetrycolors[i]
                    
                        # plot the basis elements
                        curve = PointsItem(x[boolarr], y[boolarr], size, alpha, color)
                        curve.setZValue(5)
                        graphicsview_plot[idx].addItem(curve)
            
                # --- update the graphics view ---
                graphicsview_plot[idx].repaint()
        

            
            
            
            
            
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
            self.buffer_basis = [{},{},{}]
            self.buffer_energies = [{},{},{}]
            self.buffer_positions = [{},{},{}]
            self.buffer_boolarr = [{},{},{}]
            self.buffer_basis_potential = {}
            self.buffer_energies_potential = {}
            self.buffer_positions_potential = {}
            
            self.buffer_energiesMap = [{},{},{}]
            self.buffer_positionsMap = [{},{},{}]
            self.buffer_overlapMap = [{},{},{}]
            self.buffer_energiesMap_potential = {}
            self.buffer_positionsMap_potential = {}
            self.buffer_overlapMap_potential = {}
            
            
            self.lines_buffer_minIdx = {}
            self.colormap_buffer_minIdx_potential = 0
            self.colormap_buffer_minIdx_field = [0]*3
            self.lines_buffer_minIdx_field = [0]*3
        
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
            
            # Toggle antialiasing # HACK
            if self.ui.checkbox_plot_antialiasing.isChecked():
                for plotarea in [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot, self.ui.graphicsview_potential_plot]:
                    plotarea.setAntialiasing(False)
                    plotarea.setAntialiasing(True)

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
    def toggleSamebasis(self):
        checked = self.ui.checkbox_system_samebasis.isChecked()
        if checked and self.ui.tabwidget_plotter.count() == 3:
            self.ui.tabwidget_plotter.removeTab(1)
            self.ui.tabwidget_plotter.setTabText(0, "Field map of atom 1 and 2")
        elif not checked and self.ui.tabwidget_plotter.count() == 2:
            self.ui.tabwidget_plotter.insertTab(1, self.tab_field2, "Field map of atom 2")
            self.ui.tabwidget_plotter.setTabText(0, "Field map of atom 1")
    
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
            if self.ui.checkbox_system_samebasis.isEnabled():
                self.ui.checkbox_system_samebasis.setEnabled(False)
                self.samebasis_state = self.ui.checkbox_system_samebasis.checkState()
                self.ui.checkbox_system_samebasis.setCheckState(QtCore.Qt.Unchecked) # TODO !!!!!!!!!!!
        else:
            if not self.ui.checkbox_system_samebasis.isEnabled():
                self.ui.checkbox_system_samebasis.setEnabled(True)
                self.ui.checkbox_system_samebasis.setCheckState(self.samebasis_state)
    
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
            
                if n.value()-1 < l.value():
                    n_err |= True
                    l_err |= True
            
                if abs(l.value() - j.value()) != 0.5:
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
    
    @QtCore.pyqtSlot()   
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
                    
                with open(self.path_view_last, 'w') as f:
                    params = dict()
                    params["config"] = self.ui.tabwidget_config.currentIndex()
                    params["plotter"] = self.ui.tabwidget_plotter.currentIndex()
                    params["system"] = self.ui.toolbox_system.currentIndex()
                    if self.filepath != self.userpath: params["filepath"] = self.filepath
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
                
                
                # Get limits # TODO !!!!!!!! zusammen bringen mit set limits
                self.minE_field1 = self.plotdict["minE_field1"]
                if self.minE_field1 is not None: self.minE_field1 = self.minE_field1.magnitude
                    
                self.minE_field2 = self.plotdict["minE_field2"]
                if self.minE_field2 is not None: self.minE_field2 = self.minE_field2.magnitude
                
                self.minE_potential = self.plotdict["minE_potential"]
                if self.minE_potential is not None: self.minE_potential = self.minE_potential.magnitude
                
                self.maxE_field1 = self.plotdict["maxE_field1"]
                if self.maxE_field1 is not None: self.maxE_field1 = self.maxE_field1.magnitude
                    
                self.maxE_field2 = self.plotdict["maxE_field2"]
                if self.maxE_field2 is not None: self.maxE_field2 = self.maxE_field2.magnitude
                
                self.maxE_potential = self.plotdict["maxE_potential"]
                if self.maxE_potential is not None: self.maxE_potential = self.maxE_potential.magnitude
                
            
                # Clear plots and set them up
                self.constDistance = self.getConstDistance()
                self.constEField = self.getConstEField()
                self.constBField = self.getConstBField()
                #self.sameSpecies = self.getSameSpecies()
            
                if self.senderbutton in [self.ui.pushbutton_field1_calc, self.ui.pushbutton_potential_calc] or self.samebasis:
                    self.storage_data[0] = []
                    self.storage_states[0] = None
                    self.storage_configuration[0] = [self.systemdict.paramsInOriginalunits(),self.plotdict.paramsInOriginalunits()]
                    self.storage_configuration[0][1]["gradientwidget"] = self.ui.gradientwidget_plot_gradient.saveState()
                    
                    # clear plot (with a hack)
                    autorangestate = self.ui.graphicsview_field1_plot.getViewBox().getState()["autoRange"]
                    if autorangestate[0]:
                        self.ui.graphicsview_field1_plot.disableAutoRange(axis=self.ui.graphicsview_field1_plot.getViewBox().XAxis)
                    if autorangestate[1]:
                        self.ui.graphicsview_field1_plot.disableAutoRange(axis=self.ui.graphicsview_field1_plot.getViewBox().YAxis)
                    self.ui.graphicsview_field1_plot.clear()
                    if autorangestate[0]:
                        self.ui.graphicsview_field1_plot.enableAutoRange(axis=self.ui.graphicsview_field1_plot.getViewBox().XAxis)
                    if autorangestate[1]:
                        self.ui.graphicsview_field1_plot.enableAutoRange(axis=self.ui.graphicsview_field1_plot.getViewBox().YAxis)
                        
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
                    self.leftSmallerRight = posMin < posMax
                    
                    if self.ui.checkbox_plot_autorange.isChecked():
                        self.ui.graphicsview_field1_plot.enableAutoRange()
                    else:
                        if self.ui.graphicsview_field1_plot.getViewBox().getState()["autoRange"][0]:
                            self.ui.graphicsview_field1_plot.setXRange(posMin, posMax)
                        if self.ui.graphicsview_field1_plot.getViewBox().getState()["autoRange"][1] and self.minE_field1 is not None and self.maxE_field1 is not None:
                            self.ui.graphicsview_field1_plot.setYRange(self.minE_field1, self.maxE_field1)
  
                if self.senderbutton in [self.ui.pushbutton_field2_calc, self.ui.pushbutton_potential_calc] or self.samebasis:
                    self.storage_data[1] = []
                    self.storage_states[1] = None
                    self.storage_configuration[1] = [self.systemdict.paramsInOriginalunits(),self.plotdict.paramsInOriginalunits()]
                    self.storage_configuration[1][1]["gradientwidget"] = self.ui.gradientwidget_plot_gradient.saveState()
                    
                    # clear plot (with a hack)
                    autorangestate = self.ui.graphicsview_field2_plot.getViewBox().getState()["autoRange"]
                    if autorangestate[0]:
                        self.ui.graphicsview_field2_plot.disableAutoRange(axis=self.ui.graphicsview_field2_plot.getViewBox().XAxis)
                    if autorangestate[1]:
                        self.ui.graphicsview_field2_plot.disableAutoRange(axis=self.ui.graphicsview_field2_plot.getViewBox().YAxis)
                    self.ui.graphicsview_field2_plot.clear()
                    if autorangestate[0]:
                        self.ui.graphicsview_field2_plot.enableAutoRange(axis=self.ui.graphicsview_field2_plot.getViewBox().XAxis)
                    if autorangestate[1]:
                        self.ui.graphicsview_field2_plot.enableAutoRange(axis=self.ui.graphicsview_field2_plot.getViewBox().YAxis)
                        
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
                    self.leftSmallerRight = posMin < posMax
                    
                    if self.ui.checkbox_plot_autorange.isChecked():
                        self.ui.graphicsview_field2_plot.enableAutoRange()
                    else:
                        if self.ui.graphicsview_field2_plot.getViewBox().getState()["autoRange"][0]:
                            self.ui.graphicsview_field2_plot.setXRange(posMin, posMax)
                        if self.ui.graphicsview_field2_plot.getViewBox().getState()["autoRange"][1] and self.minE_field2 is not None and self.maxE_field2 is not None:
                            self.ui.graphicsview_field2_plot.setYRange(self.minE_field2, self.maxE_field2)
                    
                if self.senderbutton == self.ui.pushbutton_potential_calc:
                    self.storage_data[2] = []
                    self.storage_states[2] = None
                    self.storage_configuration[2] = [self.systemdict.paramsInOriginalunits(),self.plotdict.paramsInOriginalunits()]
                    self.storage_configuration[2][1]["gradientwidget"] = self.ui.gradientwidget_plot_gradient.saveState()
                    
                    # clear plot (with a hack)
                    autorangestate = self.ui.graphicsview_potential_plot.getViewBox().getState()["autoRange"]
                    if autorangestate[0]:
                        self.ui.graphicsview_potential_plot.disableAutoRange(axis=self.ui.graphicsview_potential_plot.getViewBox().XAxis)
                    if autorangestate[1]:
                        self.ui.graphicsview_potential_plot.disableAutoRange(axis=self.ui.graphicsview_potential_plot.getViewBox().YAxis)
                    self.ui.graphicsview_potential_plot.clear()
                    if autorangestate[0]:
                        self.ui.graphicsview_potential_plot.enableAutoRange(axis=self.ui.graphicsview_potential_plot.getViewBox().XAxis)
                    if autorangestate[1]:
                        self.ui.graphicsview_potential_plot.enableAutoRange(axis=self.ui.graphicsview_potential_plot.getViewBox().YAxis)
                        
                    self.ui.graphicsview_potential_plot.setLabel('left', 'Energy ('+str(Units.energy)+')')
                    if self.constDistance and not self.constEField:
                        self.ui.graphicsview_potential_plot.setLabel('bottom', 'Electric field ('+str(Units.efield)+')')
                        self.converter_x_potential = Converter.fromAU(1,Units.efield).magnitude
                        posMin_potential = self.get1DPosition(self.systemdict['minEx'].magnitude,self.systemdict['minEy'].magnitude,self.systemdict['minEz'].magnitude)
                        posMax_potential = self.get1DPosition(self.systemdict['maxEx'].magnitude,self.systemdict['maxEy'].magnitude,self.systemdict['maxEz'].magnitude)
                    elif self.constDistance and not self.constBField:
                        self.ui.graphicsview_potential_plot.setLabel('bottom', 'Magnetic field ('+str(Units.bfield)+')')
                        self.converter_x_potential = Converter.fromAU(1,Units.bfield).magnitude
                        posMin_potential = self.get1DPosition(self.systemdict['minBx'].magnitude,self.systemdict['minBy'].magnitude,self.systemdict['minBz'].magnitude)
                        posMax_potential = self.get1DPosition(self.systemdict['maxBx'].magnitude,self.systemdict['maxBy'].magnitude,self.systemdict['maxBz'].magnitude)
                    else:
                        self.ui.graphicsview_potential_plot.setLabel('bottom', 'Interatomic distance ('+str(Units.length)+')')
                        self.converter_x_potential = Converter.fromAU(1,Units.length).magnitude
                        posMin_potential = self.systemdict['minR'].magnitude
                        posMax_potential = self.systemdict['maxR'].magnitude
                    self.leftSmallerRight_potential = posMin_potential < posMax_potential
                    
                    if self.ui.checkbox_plot_autorange.isChecked():
                        self.ui.graphicsview_potential_plot.enableAutoRange()
                    else:
                        if self.ui.graphicsview_potential_plot.getViewBox().getState()["autoRange"][0]:
                            self.ui.graphicsview_potential_plot.setXRange(posMin_potential, posMax_potential)
                        if self.ui.graphicsview_potential_plot.getViewBox().getState()["autoRange"][1] and self.minE_potential is not None and self.maxE_potential is not None:
                            self.ui.graphicsview_potential_plot.setYRange(self.minE_potential, self.maxE_potential)
                        
                self.converter_y = Converter.fromAU(1,Units.energy).magnitude
                
                
                self.steps = self.systemdict['steps']
                
                
                
                n1 = self.systemdict['n1']
                l1 = self.systemdict['l1']
                j1 = self.systemdict['j1']
                m1 = self.systemdict['m1']
                n2 = self.systemdict['n2']
                l2 = self.systemdict['l2']
                j2 = self.systemdict['j2']
                m2 = self.systemdict['m2']
                self.unperturbedstate = np.array([n1,l1,j1,m1,n2,l2,j2,m2])
                
                if self.ui.groupbox_plot_overlap.isChecked():
                    if self.plotdict["overlapUnperturbed"]:
                        self.overlapstate = self.unperturbedstate
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
        
        if idx == 0 and self.ui.checkbox_system_samebasis.isChecked(): description = "field map of atom 1 and 2"
        else: description = ["field map of atom 1", "field map of atom 2", "pair potential"][idx]
        
        filename,_ = QtGui.QFileDialog.getSaveFileName(self, \
            "Save {}".format(description), path, "zip (*.zip)")
        
        if not filename:
            return
        
        self.resultfile = filename
        self.filepath = os.path.dirname(filename)
        
        # open zip file
        ziparchive = zipfile.ZipFile(filename, 'w', compression = zipfile.ZIP_STORED) # zipfile.ZIP_DEFLATED
        
        try:
            # save plot
            plotitem = [self.ui.graphicsview_field1_plot, self.ui.graphicsview_field2_plot, self.ui.graphicsview_potential_plot][idx].getPlotItem()
            exporter = pg.exporters.ImageExporter(plotitem)
            exporter.parameters()['width'] = 2000
            exporter.parameters()['height'] = 2000
            exporter.parameters()['antialias'] = True
            image = exporter.export(toBytes=True)
            
            buffer = QtCore.QBuffer()
            buffer.open(QtCore.QIODevice.WriteOnly)
            image = image.scaled(1500, 1500, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
            image.save(buffer, "PNG")
            ziparchive.writestr('plot.png', buffer.data())
            
            # save configuration
            filelike = StringIO()
            json.dump(self.storage_configuration[idx][0], filelike, indent=4, sort_keys=True)
            ziparchive.writestr('settings.sconf', filelike.getvalue())
            
            filelike = StringIO()
            json.dump(self.storage_configuration[idx][1], filelike, indent=4, sort_keys=True)
            ziparchive.writestr('settings.pconf', filelike.getvalue())
            
            # create data dictionary
            data = {}
            
            data['numStates'] = 0
            data['numSteps'] = 0
            data['numEigenpairs'] = 0
            data['states'] = []
            data['eigenvectors'] = []
            data['eigenvalues'] = []
            data['bfields'] = []
            data['efields'] = []
            if idx == 2: data['distances'] = []
            
            if idx == 0 or idx == 1: data['states_description'] = 'state(idxState, {idxState, n, l, j, m})'
            elif idx == 2: data['states_description'] = 'state(idxState, {idxState, n1, l1, j1, m1, n2, l2, j2, m2})'
            data['eigenvectors_description'] = 'eigenvectorcoordinate(0, idxStep, idxState, idxEigenpair)'
            data['eigenvalues_description'] = 'eigenvalue(idxStep, idxEigenpair)'
            data['bfields_description'] = 'bfieldcoordinate(idxStep, {Bx, By, Bz})'
            data['efields_description'] = 'efieldcoordinate(idxStep, {Ex, Ey, Ez})'
            if idx == 2: data['distances_description'] = 'distance(0, idxStep)'
            
            # save states
            if self.storage_states[idx] is not None:
                data['states'] = self.storage_states[idx] # nState, i-n1-l1-j1-m1-n2-l2-j2-m2 # nState, i-n-l-j-m
                data['numStates'] = len(data['states'])
            
            # save data
            self.converter_bfield = Converter.fromAU(1,Units.bfield).magnitude # TODO Variable an anderer Stelle anlegen
            self.converter_efield = Converter.fromAU(1,Units.efield).magnitude # TODO Variable an anderer Stelle anlegen
            self.converter_length = Converter.fromAU(1,Units.length).magnitude # TODO Variable an anderer Stelle anlegen
            
            filestep_last = None
            
            for filestep, blocknumber, filename in sorted(self.storage_data[idx],key=itemgetter(0,1)):
                eigensystem = Eigensystem(filename)
                energies = eigensystem.energies*self.converter_y # nBasis
                basis = eigensystem.basis  # nState, nBasis (stored in Compressed Sparse Column format, CSC)
                    
                if filestep != filestep_last: # new step
                    data['bfields'].append(np.array([float(eigensystem.params["Bx"]),float(eigensystem.params["By"]),float(eigensystem.params["Bz"])])*self.converter_efield)
                    data['efields'].append(np.array([float(eigensystem.params["Ex"]),float(eigensystem.params["Ey"]),float(eigensystem.params["Ez"])])*self.converter_bfield)
                    if idx == 2: data['distances'].append(float(eigensystem.params["R"])*self.converter_length)
                    data['eigenvalues'].append(energies)
                    data['eigenvectors'].append(basis)
                else: # new block
                    data['eigenvalues'][-1] = np.append(data['eigenvalues'][-1],energies)
                    csc_happend(data['eigenvectors'][-1],basis)
                    
                filestep_last = filestep  
            
            if len(data['eigenvalues']) > 0:
                data['numSteps'] = len(data['eigenvalues'])
                data['numEigenpairs'] = len(data['eigenvalues'][0])
            
            filelike=BytesIO()
            io.savemat(filelike,data,do_compression=False,format='5',oned_as='row')
            ziparchive.writestr('data.mat', filelike.getvalue())
                            
        finally:
            # close zip file
            ziparchive.close()
        
        
        """ symmetry = eigensystem.params["symmetry"] # TODO ausgelesene Symmetrie auch beim Plotten verwenden
                
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
        filename,_ = QtGui.QFileDialog.getSaveFileName(self, \
            "Save system configuration",path, "sconf (*.sconf)")
        
        if filename:
            with open(filename, 'w') as f:
                params = self.systemdict.paramsInOriginalunits()
                json.dump(params, f, indent=4, sort_keys=True)
            self.systemfile = filename
            self.filepath = os.path.dirname(filename)
    
    @QtCore.pyqtSlot()
    def savePlotConf(self):
        path = self.plotfile if self.plotfile is not None else self.filepath
        filename,_ = QtGui.QFileDialog.getSaveFileName(self, \
            "Save plot configuration",path, "pconf (*.pconf)")
        
        if filename:
            with open(filename, 'w') as f:
                params = self.plotdict.paramsInOriginalunits()
                params["gradientwidget"] = self.ui.gradientwidget_plot_gradient.saveState()
                json.dump(params, f, indent=4, sort_keys=True)
            self.plotfile = filename
            self.filepath = os.path.dirname(filename)
    
    @QtCore.pyqtSlot()
    def openSystemConf(self):
        filename,_ = QtGui.QFileDialog.getOpenFileName(self, \
            "Open system configuration",self.filepath, "sconf (*.sconf)")
        
        if filename:
            with open(filename, 'r') as f:
                params = json.load(f)
                self.systemdict.load(params)
            self.systemfile = filename
            self.filepath = os.path.dirname(filename)
    
    @QtCore.pyqtSlot()
    def openPlotConf(self):
        filename,_ = QtGui.QFileDialog.getOpenFileName(self, \
            "Open plot configuration",self.filepath, "pconf (*.pconf)")
        
        if not (filename == ""):
            with open(filename, 'r') as f:
                params = json.load(f)
                self.ui.gradientwidget_plot_gradient.restoreState(params["gradientwidget"])
                del params["gradientwidget"]
                self.plotdict.load(params)
            self.plotfile = filename
            self.filepath = os.path.dirname(filename)
            
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
        
        with open(self.path_view_last, 'w') as f:
            params = dict()
            params["config"] = self.ui.tabwidget_config.currentIndex()
            params["plotter"] = self.ui.tabwidget_plotter.currentIndex()
            params["system"] = self.ui.toolbox_system.currentIndex()
            if self.filepath != self.userpath: params["filepath"] = self.filepath
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
