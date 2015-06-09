#!/usr/bin/env python

from version import *
import system_settings as ss

import wx
import wx.lib.scrolledpanel
import os
import json

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

class ComboBox(wx.ComboBox):
    def SetValue(self,val):
        self.SetStringSelection(val)

    def GetValue(self):
        return self.GetCurrentSelection()

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
        gSizer = wx.FlexGridSizer(rows=len(ss.par), cols=2)
        gSizer.AddGrowableCol(0, 0)
        gSizer.AddGrowableCol(1, 1)

        #gSizer.Add(ComboBox(self, wx.ID_ANY, 'ScaLAPACK', choices=['SLEPc_KrylovSchur', 'SLEPc_Lapack', 'ScaLAPACK'], style=wx.CB_READONLY), 0, wx.EXPAND)
        
        for key in ss.par:
            self.__addLabel(gSizer, ss.par[key]['name'])
            
            ss.par[key]['wx'] = wx.TextCtrl(self, wx.ID_ANY, str(ss.par[key]['value']))
            gSizer.Add(ss.par[key]['wx'], 0, wx.EXPAND)

        vSizer.Add(gSizer, 0, wx.ALL|wx.EXPAND, 5)
        self.SetSizer(vSizer)
        vSizer.Fit(self)

    def __addLabel(self, gSizer, label):
        label = wx.StaticText(self, wx.ID_ANY, label)
        hSizer = wx.BoxSizer(wx.HORIZONTAL)
        hSizer.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
        gSizer.Add(hSizer, 0, wx.ALIGN_RIGHT)
        

class PlotPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        # Initialize matplotlib
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.updatePlot()
        
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

    def updatePlot(self):
        self.axes.plot([1,3,2])


class Notebook(wx.Notebook):
    def __init__(self, parent):
        wx.Notebook.__init__(self, parent, wx.ID_ANY, style=wx.BK_DEFAULT)
        self.AddPage(InputPanel(self), "Parameters")
        self.AddPage(PlotPanel(self), "Plot")


class MyForm(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, title=APP_NAME,size=(600,600))

        # StatusBar
        self.CreateStatusBar()

        # MenuBar
        self.__createMenu()

        # Panel
        panel = wx.Panel(self)
 
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(Notebook(panel), 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)
        self.Layout()

        self.Bind(wx.EVT_CLOSE, self.__OnExit)
 
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

    def __OnExit(self, e):
        print "Close window called"
        self.Destroy()

    def __OnOpen(self, e):
        self.dirname = ''
        dlg = wx.FileDialog(self, "Load file", self.dirname, "", "*.json", wx.OPEN)
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
        dlg = wx.FileDialog(self, "Save file", self.dirname, "", "*.json", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            f = open(os.path.join(self.dirname, self.filename), 'w')
            parameters = {}
            for key in ss.par:
                try:
                    parameters[key] = ss.par[key]['type'](ss.par[key]['wx'].GetValue())
                except:
                    print "Type Error!"
            json.dump(parameters, f, indent=4)
            f.close()
        dlg.Destroy()


def main():
    app = wx.App(False)
    frame = MyForm()
    frame.Show(True)
    app.MainLoop()
