#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
from matplotlib import rcParams
font = {'family' : 'Dejavu Sans',
        'weight' : 'normal',
        'size'   : 30}
rc('font', **font)
rcParams['lines.linewidth'] = 4
rcParams['lines.markersize'] = 18
rcParams['markers.fillstyle'] = 'none'
rcParams.update({'figure.autolayout': True})


def plotLogLogData (data_x, data_y, ref_data_y, myPlotDict):
    
    fig, ax = plt.subplots()
    for k in range(0,len(data_y)):
        N = data_y[k].shape[0]
        ax.loglog(data_x[k][0:N], data_y[k], myPlotDict['data_markers'][k], label=myPlotDict['data_labels'][k])
    
    for k in range(0,len(ref_data_y)):
        N = ref_data_y[k].shape[0]
        ax.loglog(data_x[k][0:N], ref_data_y[k], myPlotDict['ref_data_markers'][k], label=myPlotDict['ref_data_labels'][k])
    
    legend = ax.legend(loc=myPlotDict['legend_loc'], shadow=False)
    plt.xlabel(myPlotDict['xlabel'])
    plt.ylabel(myPlotDict['ylabel'])
    plt.xlim(myPlotDict['xlim'])
    # plt.ylim(myPlotDict['ylim'])
    plt.grid(which='major', linestyle='-', linewidth=3)
    #plt.grid(which='minor', linestyle='-.', linewidth=0.1)
    plt.title(myPlotDict['title'])
    
    fig = plt.gcf()
    fig.set_size_inches(18, 14)
    if myPlotDict['show_plot'] and not myPlotDict['save_plot']:
        plt.show()
    elif myPlotDict['save_plot']:
        fig.savefig(myPlotDict['out_filename'], format='pdf', dpi=1000)
    plt.clf()
    

def plotData (data_x, data_y, myPlotDict):
    
    fig, ax = plt.subplots()
    for k in range(0,len(data_y)):
        N = data_y[k].shape[0]
        ax.plot(data_x[k][0:N], data_y[k], myPlotDict['data_markers'][k], label=myPlotDict['data_labels'][k])
    
    plt.xlabel(myPlotDict['xlabel'])
    plt.ylabel(myPlotDict['ylabel'])
    plt.xlim(myPlotDict['xlim'])
    #plt.ylim(myPlotDict['ylim'])
    plt.grid()
    plt.title(myPlotDict['title'])
    
    fig = plt.gcf()
    fig.set_size_inches(18, 14)
    if myPlotDict['show_plot'] and not myPlotDict['save_plot']:
        plt.show()
    elif myPlotDict['save_plot']:
        fig.savefig(myPlotDict['out_filename'], format='pdf', dpi=1000)
    plt.clf()

## End of file
