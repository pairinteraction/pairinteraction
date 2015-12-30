from version import *
import system_settings as ss

import wx
import wx.lib.scrolledpanel
import os
import json

import subprocess, signal, threading, multiprocessing
import numpy as np
from scipy import sparse

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

# --- methods to read energies and basis from files ---

typeIds = {1008: 'int8', 1016 : 'int16', 1032 : 'int32', 1064 : 'int64', 1108 : 'uint8', 1116 : 'uint16', 1132 : 'uint32', \
    1164 : 'int64', 2032 : 'float32', 2064 : 'float64'}
type_t = 'uint16'
csr_not_csc = 0x01; # xxx0: csc, xxx1: csr
complex_not_real = 0x02; # xx0x: real, xx1x: complex

def readNumber(f, sz = 1):
    datatype = typeIds[np.fromfile(f, dtype=np.dtype(type_t), count=1)[0]]
    return np.squeeze(np.fromfile(f, dtype=np.dtype(datatype), count=sz))

def readVector(f):
    size = readNumber(f)
    return readNumber(f, size)

def readMatrix(f):
    flags = readNumber(f)
    rows = readNumber(f)
    cols = readNumber(f)
    if flags & complex_not_real: data = readVector(f) + readVector(f)*1j
    else: data = readVector(f)
    indices = readVector(f)
    indptr = np.append(readVector(f),len(data))
    if flags & csr_not_csc: return sparse.csr_matrix((data, indices, indptr), shape=(rows, cols))
    else: return sparse.csc_matrix((data, indices, indptr), shape=(rows, cols))

def load(filename):
    with open(filename,'rb') as f:
        energies = readMatrix(f).diagonal()
        basis = readMatrix(f)
        return energies, basis

def load2(filename):
    with open(filename,'rb') as f:
        energies = readMatrix(f)
        return energies
        
# --- definition of the gui ---

# IDs
ID_STARTABORT = wx.NewId()
EVT_RESULT_ID = wx.NewId()
EVT_FINISH_ID = wx.NewId()

# Define notification event for thread completion
def EVT_RESULT(win, func):
    win.Connect(-1, -1, EVT_RESULT_ID, func)

def EVT_FINISH(win, func):
    win.Connect(-1, -1, EVT_FINISH_ID, func)

class ResultEvent(wx.PyEvent):
    def __init__(self, data):
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_RESULT_ID)
        self.data = data

class FinishEvent(wx.PyEvent):
    def __init__(self):
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_FINISH_ID)

# Thread class
class WorkerThread(threading.Thread):
    def __init__(self, notify_window, stdout):
        threading.Thread.__init__(self)
        self._notify_window = notify_window
        self._stdout = stdout
        self.start()

    def run(self):
        for line in iter(self._stdout.readline, b""):
            #print (line.decode('utf-8'), end="")
            if line[:5] == b">>OUT":
                wx.PostEvent(self._notify_window, ResultEvent(line[6:-1].decode('utf-8')))
            if line[:5] == b">>END" or not line:
                break
        
        wx.PostEvent(self._notify_window, FinishEvent())

# Combo box
class ComboBox(wx.ComboBox):
    def SetValue(self,val):
        self.SetStringSelection(val)

    def GetValue(self):
        return self.GetCurrentSelection()

# Input panel
class InputPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent):
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent)
        self.SetupScrolling()

        # Sizer
        vSizer = wx.BoxSizer(wx.VERTICAL)
        hSizer = wx.BoxSizer(wx.HORIZONTAL)
        hSizer.Add(wx.StaticText(self, wx.ID_ANY, APP_NAME), 0, wx.ALL , 5)
        vSizer.Add(hSizer, 0, wx.CENTER)
        vSizer.Add(wx.StaticLine(self,), 0, wx.ALL|wx.EXPAND, 5)

        # Maybe replace with GridBagSizer
        gSizer = wx.FlexGridSizer(rows=len(ss.par), cols=2, gap=wx.Size(0,0))
        gSizer.AddGrowableCol(0, 0)
        gSizer.AddGrowableCol(1, 1)

        #gSizer.Add(ComboBox(self, wx.ID_ANY, 'ScaLAPACK', choices=['SLEPc_KrylovSchur', 'SLEPc_Lapack', 'ScaLAPACK'], style=wx.CB_READONLY), 0, wx.EXPAND)
        
        for key in ss.par:
            self.__addLabel(gSizer, ss.par[key]['name'])
            
            ss.par[key]['wx'] = wx.TextCtrl(self, wx.ID_ANY, str(ss.par[key]['value']))
            gSizer.Add(ss.par[key]['wx'], 0, wx.EXPAND)
            
        # button
        self.buttonStartAbort = wx.Button(self, ID_STARTABORT, 'Start calculation',pos=(0,0)) # TODO vernuenftig ausrichten

        vSizer.Add(gSizer, 0, wx.ALL|wx.EXPAND, 5)
        self.SetSizer(vSizer)
        vSizer.Fit(self)

    def __addLabel(self, gSizer, label):
        label = wx.StaticText(self, wx.ID_ANY, label)
        hSizer = wx.BoxSizer(wx.HORIZONTAL)
        hSizer.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
        gSizer.Add(hSizer, 0, wx.ALIGN_RIGHT)
        
# Plot panel
class PlotPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        # Initialize matplotlib
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.initPlot()
        
        # Bind matplotlib to wx
        self.canvas = FigureCanvasWxAgg(self, wx.ID_ANY, self.figure)
        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.toolbar.Realize()

        # Arrange binded objects in window
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Fit()

    def initPlot(self):        
        self.axes.set_xlim(-0.005,0.05)
        self.axes.set_ylim(2.0,4.8)
    
    def append(self, x, y):
    
        au2MHz = 6579683.920729
        au2Vpercm = 5.14220652e11/100
        
        x *= au2Vpercm
        y *= au2MHz
        y = y[(y <= 4.8) & (y >= 2.0)]

        self.axes.plot(x,y.reshape(1,-1),'kx',ms=2)
        self.figure.canvas.draw()

# Notebook containing the panels
class Notebook(wx.Notebook):
    def __init__(self, parent):
        wx.Notebook.__init__(self, parent, wx.ID_ANY, style=wx.BK_DEFAULT)
        self.input = InputPanel(self)
        self.plot = PlotPanel(self)
        self.AddPage(self.input, "Parameters")
        self.AddPage(self.plot, "Plot")

# GUI class
class MyForm(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title=APP_NAME,size=(600,600))

        # StatusBar
        self.CreateStatusBar()

        # MenuBar
        self.__createMenu()

        # Panel
        panel = wx.Panel(self)
 
        self.notebook = Notebook(panel)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)
        self.Layout()

        self.Bind(wx.EVT_CLOSE, self.__OnExit)
        
        self.Bind(wx.EVT_BUTTON, self.__OnStartAbort, id=ID_STARTABORT)
        
        # Set up event handler for any worker thread results
        EVT_RESULT(self,self.__OnResult)
        EVT_FINISH(self,self.__OnFinish)
        
        # Indicate that there is no c++ process running yet
        self.p = None
        
        # Indicate that there is no thread running yet
        self.worker = None
        
        self.Show()
        self.SetMinSize((400,400))


    def __createMenu(self):
        # Create and populate menu
        filemenu = wx.Menu()
        menuOpen = filemenu.Append(wx.ID_OPEN, "&Load", "Load file")
        menuSave = filemenu.Append(wx.ID_SAVE, "&Save", "Save file")
        menuExit = filemenu.Append(wx.ID_EXIT, "&Exit", "Exit program")

        # MenuBar
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File")
        self.SetMenuBar(menuBar)

        # Bind menu events
        self.Bind(wx.EVT_MENU, self.__OnOpen, menuOpen)
        self.Bind(wx.EVT_MENU, self.__OnSave, menuSave)
        self.Bind(wx.EVT_MENU, self.__OnExit, menuExit)
    
    def __abort(self, e):
        self.p.kill()
        print("Calculation aborted")

    def __start(self, e):
        # Save configuration to json file
        f = open("params.json", 'w')
        parameters = {}
        for key in ss.par:
            try:
                parameters[key] = ss.par[key]['type'](ss.par[key]['wx'].GetValue())
            except:
                print("Type Error!")
        json.dump(parameters, f, indent=4)
        f.close()
        
        # Start c++ process
        numprocessors   = multiprocessing.cpu_count()
        path_workingdir = "../build/"
        path_cpp = "../build/pairinteraction"
        path_config = "../gui/params.json"
        
        self.p = subprocess.Popen(["mpiexec","-n","%d"%numprocessors,path_cpp,"-c",path_config],
            stdout=subprocess.PIPE, cwd=path_workingdir)
        
        # Start thread that collects the output
        self.worker = WorkerThread(self, self.p.stdout)
        
        print("Calculation started")

    def __OnExit(self, e):
        print("Close window called")
        if self.p is not None:
            self.__abort(e)
        self.Destroy()
    
    def __OnStartAbort(self, e):
        if self.p is not None:
            self.__abort(e)
        else:
            self.__start(e)
            self.notebook.input.buttonStartAbort.SetLabel("Abort calculation")

    def __OnResult(self, event):
        print(event.data)
    
        path_json = event.data + ".json"
        path_mat = event.data + ".mat"
        
        with open(path_json, 'r') as f:
            params = json.load(f)
        
            vecB = np.array([params["Bx"],params["By"],params["Bz"]]).astype(float)
            vecE = np.array([params["Ex"],params["Ey"],params["Ez"]]).astype(float)
            normB = np.linalg.norm(vecB)
            normE = np.linalg.norm(vecE)
            
            energies = load(path_mat)[0].real
            
            self.notebook.plot.append(normE,energies) # TODO
                
    def __OnFinish(self, event):
        self.p = None
        self.notebook.input.buttonStartAbort.SetLabel("Start calculation")

    def __OnOpen(self, e):
        self.dirname = ''
        dlg = wx.FileDialog(self, "Load file", self.dirname, "", "*.json", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            f = open(os.path.join(self.dirname, self.filename), 'r')
            parameters = json.load(f)
            for key in parameters:
                ss.par[key]['value'] = parameters[key]
                ss.par[key]['wx'].SetValue(str(parameters[key]))
            f.close()
        dlg.Destroy()

    def __OnSave(self, e):
        self.dirname = ''
        dlg = wx.FileDialog(self, "Save file", self.dirname, "", "*.json", wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            f = open(os.path.join(self.dirname, self.filename), 'w')
            parameters = {}
            for key in ss.par:
                try:
                    parameters[key] = ss.par[key]['type'](ss.par[key]['wx'].GetValue())
                except:
                    print("Type Error!")
            json.dump(parameters, f, indent=4)
            f.close()
        dlg.Destroy()


def main():
    app = wx.App(False)
    frame = MyForm()
    frame.Show(True)
    app.MainLoop()
