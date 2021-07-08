#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import numpy as np
import common


def load_json(filename):
    with open(filename, "r") as f:
        return json.load(f)

def load_data(filename):
    j = load_json(filename)
    return np.array(j["time"]), np.array(j["pressure"])


def pressure_signal(dir_sol, dir_fig, lx, show_plot, save_plot):

    fn = dir_sol+'lx'+str(lx)+'/pressure_signal_lx'+str(lx)+'.json'
    time, pressure = load_data(fn)
    
    N = pressure.size
    signal = np.zeros(N)
    for k in range(1, N):
        h = time[k] - time[k - 1]
        signal[k] = signal[k - 1] + 0.5 * (pressure[k] + pressure[k - 1]) * h
    
    myPlotDict = {}
    myPlotDict['show_plot'] = show_plot
    myPlotDict['save_plot'] = save_plot
    
    myPlotDict['legend_loc'] = 'lower right'
    myPlotDict['data_markers'] = ['b.-']
    
    myPlotDict['title'] = ''
    myPlotDict['xlabel'] = r'time $t$'
    myPlotDict['ylabel'] = r'$u_{C}(t)$'
    myPlotDict['data_labels'] = ['']
    myPlotDict['out_filename'] =dir_fig+'scattering_signal_lx'+str(lx)+'.pdf'
    myPlotDict['xlim'] = [-0.1, 1.1]
    myPlotDict['ylim'] = []
    common.plotData([time], [signal], myPlotDict)


if __name__ == "__main__":

    show_plot = True
    save_plot = True
    
    dir_sol = '../../output/transmission_test1/'
    dir_fig = '../../figures/transmission_test1/'
    
    lx = 5
    pressure_signal(dir_sol, dir_fig, lx, show_plot, save_plot)
    
    
# End of file
