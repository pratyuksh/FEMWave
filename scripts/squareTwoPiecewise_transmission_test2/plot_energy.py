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
    return np.array(j["time"]), np.array(j["energy"])


def energy_evolution(dir_sol, dir_fig, lx, show_plot, save_plot):

    fn1 = dir_sol+'qu/lx'+str(lx)+'/energy_evolution_lx'+str(lx)+'.json'
    time1, energy1 = load_data(fn1)
    
    fn2 = dir_sol+'br/lx'+str(lx)+'/energy_evolution_lx'+str(lx)+'.json'
    time2, energy2 = load_data(fn2)
    
    myPlotDict = {}
    myPlotDict['show_plot'] = show_plot
    myPlotDict['save_plot'] = save_plot
    
    myPlotDict['legend_loc'] = 'lower right'
    myPlotDict['data_markers'] = ['b.-', 'r.-']
    
    myPlotDict['title'] = ''
    myPlotDict['xlabel'] = r'time $t$'
    myPlotDict['ylabel'] = r'Energy $\mathcal{E}(t)$'
    myPlotDict['data_labels'] = ['qu', 'br']
    myPlotDict['out_filename'] =dir_fig+'energy_evolution_lx'+str(lx)+'.pdf'
    myPlotDict['xlim'] = [-0.02, 0.32]
    myPlotDict['ylim'] = []
    common.plotData([time1, time2], [energy1, energy2], myPlotDict)


if __name__ == "__main__":

    show_plot = True
    save_plot = False
    
    dir_sol = '../../output/squareTwoPiecewise_transmission_test2/'
    dir_fig = '../../figures/squareTwoPiecewise_transmission_test2/'
    
    lx = 5
    energy_evolution(dir_sol, dir_fig, lx, show_plot, save_plot)
    
    
# End of file
