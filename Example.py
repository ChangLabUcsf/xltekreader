
"""
Example.py

Last Edited: 12/1/2017

Lead Author[s]: Anthony Fong
Contributor[s]:

Description:


Machine I/O
Input:
Output:

User I/O
Input:
Output:


"""
########################################################################################################################

########## Libraries, Imports, & Setup ##########

# Default Libraries #
import os
import sys
import time
import pathlib
import datetime
import multiprocessing

# Downloaded Libraries #
import h5py
import numpy as np

from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style

# Imports from Local Packages #
import xltek_reader






# Classes #

class XLTEK_Plotter:
    def __init__(self, subj, reader=None):
        self.subj = subj
        if reader is None:
            self.reader = xltek_reader.XLTEK_Reader(subj, update=False)
        else:
            self.reader = reader

        self.buffer = None

        style.use('fivethirtyeight')
        self.x_plots = 5
        self.y_plots = 9
        self.fig, self.axes = plt.subplots(self.y_plots, self.x_plots)
        # self.lines = self.axes.plot([],[],'-',animated=True)
        # self.fig, self.axes = plt.subplots(5, 5)
        self.lines = []
        for i in range(0,self.y_plots):
            for j in range(0,self.x_plots):
                self.axes[i,j].set_title('Channel %d'%(i*self.x_plots+j+1),fontsize=5)
                self.axes[i,j].set_xlim(0, 20000)
                xticks = [0, 5000, 10000, 15000, 20000]
                self.axes[i, j].set_xticks(xticks)
                self.axes[i, j].set_xticklabels(["%d" % (x/1000-20) for x in xticks])
                self.axes[i,j].xaxis.set_tick_params(labelsize=5)
                self.axes[i,j].set_ylim(-1000, 1000)
                yticks = [-1000, -750, -500, -250, 0, 250, 500, 750, 1000]
                self.axes[i, j].set_yticks(yticks)
                self.axes[i, j].set_yticklabels(["%d" % y for y in yticks])
                self.axes[i, j].yaxis.set_tick_params(labelsize=5)
                self.lines += self.axes[i,j].plot([], [], '-', linewidth=0.25, animated=True)
        for i in range(0, self.y_plots-1):
            plt.setp([a.get_xticklabels() for a in self.axes[i, :]], visible=False)
        for i in range(1,self.x_plots):
            plt.setp([a.get_yticklabels() for a in self.axes[:, i]], visible=False)

    def first_frame(self):
        self.axes.set_xlim(0,3000)
        self.axes.set_ylim(-100,100)
        return self.lines,

    def spawn_collector(self):
        self.buffer = self.reader.start_current()
        return self.buffer

    def update(self, frame):
        data = self.buffer[-20:]
        ecog = np.concatenate(tuple(entry['data'] for entry in data), axis=0)
        dims = ecog.shape
        xs = range(0, dims[0])
        #self.line.set_data(xs, ecog[:, 20])
        for i in range(0,self.y_plots):
            for j in range(0,self.x_plots):
                self.lines[i*self.x_plots+j].set_data(xs, ecog[:,i*self.x_plots+j])
        return self.lines

    def animate(self):
        self.animator = animation.FuncAnimation(self.fig, self.update, interval=1000, blit=True)
        plt.show()

    def save(self):
        Writer = animation.writers['ffmpeg']
        self.writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        self.animator = animation.FuncAnimation(self.fig, self.update, interval=50, blit=True)
        self.animator.save('lines.mp4', writer=self.writer)
        plt.show()

class Lowry_Preprocess:
    def __init__(self, subj, reader=None):
        self.subj = subj
        if reader is None:
            self.reader = xltek_reader.XLTEK_Reader(subj, update=False)
        else:
            self.reader = reader

        self.buffer = None

        self.stop_events = {}

        # pyqtgraph stuff
        pg.setConfigOptions(antialias=True)
        self.traces = dict()
        self.app = QtGui.QApplication(sys.argv)
        self.win = pg.GraphicsWindow(title='Lowry Preprocessed')
        self.win.setWindowTitle('Lowry Preprocessed')
        self.win.setGeometry(5, 115, 1910, 1070)

        self.x_plots = 6
        self.y_plots = 9
        self.axes = [[] for p in range(0,self.y_plots)]
        self.plots = [[] for p in range(0,self.y_plots)]
        for i in range(0, self.y_plots):
            for j in range(0, self.x_plots):
                title = ('Channel %d'%(i * self.x_plots + j + 1))
                wf_xaxis = pg.AxisItem(orientation='bottom')
                wf_xaxis.setStyle(autoExpandTextSpace=False)

                wf_yaxis = pg.AxisItem(orientation='left')
                self.axes[i].append(self.win.addPlot(title=title, row=i, col=j, axisItems={'bottom': wf_xaxis, 'left': wf_yaxis}))
                # self.axes[i][j].setYRange(-1000, 1000, padding=0)
                self.axes[i][j].setXRange(-1000, 21000, padding=0.005)
                self.plots[i].append(self.axes[i][j].plot(pen='c', width=1))

    def spawn_collector(self):
        self.buffer = self.reader.start_current()
        return self.buffer

    def spawn_forward(self, start_entry):
        self.buffer, self.forward_process, self.stop_events['buffer'] = self.reader.start_forward(start_entry, spawn=True)
        return self.buffer

    def update(self):
        data = self.buffer[-20:]
        ecog = np.concatenate(tuple(entry['data'] for entry in data), axis=0)
        dims = ecog.shape
        xs = range(0, dims[0])
        for i in range(0, self.y_plots):
            for j in range(0, self.x_plots):
                self.plots[i][j].setData(xs,ecog[:,i*self.x_plots+j])
                wf_xlabels = [(0, data[0]['snc_start'].time().isoformat(timespec='seconds')),
                              (10000, data[10]['snc_start'].time().isoformat(timespec='seconds')),
                              (len(ecog), data[-1]['snc_end'].time().isoformat(timespec='seconds'))]
                self.axes[i][j].getAxis('bottom').setTicks([wf_xlabels])
                self.axes[i][j].showAxis('bottom')

    def for_update(self):
        data = self.buffer[:20]
        self.buffer.pop(0)
        ecog = np.concatenate(tuple(entry['data'] for entry in data), axis=0)
        dims = ecog.shape
        xs = range(0, dims[0])
        for i in range(0, self.y_plots):
            for j in range(0, self.x_plots):
                self.plots[i][j].setData(xs, ecog[:, i * self.x_plots + j])
                wf_xlabels = [(0, data[0]['snc_start'].time().isoformat(timespec='seconds')),
                              (10000, data[10]['snc_start'].time().isoformat(timespec='seconds')),
                              (len(ecog), data[-1]['snc_end'].time().isoformat(timespec='seconds'))]
                self.axes[i][j].getAxis('bottom').setTicks([wf_xlabels])
                self.axes[i][j].showAxis('bottom')

    def start(self):
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtGui.QApplication.instance().exec_()

    def for_animation(self):
        timer = QtCore.QTimer()
        timer.timeout.connect(self.for_update)
        timer.start(1000)
        self.start()

    def animation(self):
        timer = QtCore.QTimer()
        timer.timeout.connect(self.update)
        timer.start(20)
        self.start()

if __name__ == '__main__':

    thing = Lowry_Preprocess('test',[])
    time.sleep(30)