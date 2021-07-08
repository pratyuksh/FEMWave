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
    return np.array(j["h_max"]), np.array(j["L2Error"]), np.array(j["H1Error"])


def convg(dir_sol, dir_fig, show_plot, save_plot):

    fn1 = dir_sol+'convergence_lShaped_singular_interpolant_qu.json'
    fn2 = dir_sol+'convergence_lShaped_singular_interpolant_br.json'
    hMax, L2err1, H1err1 = load_data(fn1)
    hMax, L2err2, H1err2 = load_data(fn2)
    
    err1 = np.sqrt(L2err1**2 + H1err1**2)
    err2 = np.sqrt(L2err2**2 + H1err2**2)
    print('Error Qu: ', err1)
    print('Error Br: ', err2)
    
    # linear fitting
    m = hMax.size
    linfit1 = np.polyfit(np.log(hMax[m-3:m]), np.log(err1[m-3:m]), 1)
    linfit2 = np.polyfit(np.log(hMax[m-3:m]), np.log(err2[m-3:m]), 1)
    ref1 = np.exp(np.polyval(linfit1, np.log(hMax))+0.5)
    ref2 = np.exp(np.polyval(linfit2, np.log(hMax))-0.5)
    slope1 = linfit1[0]
    slope2 = linfit2[0]
    
    myPlotDict = {}
    myPlotDict['show_plot'] = show_plot
    myPlotDict['save_plot'] = save_plot
    
    # plot error vs mesh size
    myPlotDict['xlabel'] = r'Meshsize $h$ [log]'
    myPlotDict['legend_loc'] = 'upper left'
    myPlotDict['data_markers'] = ['bs-', 'ro-']
    myPlotDict['data_labels'] = ['qu', 'br']
    myPlotDict['ylim'] = []
    myPlotDict['ref_data_markers'] = ['b--', 'r-.']
    
    myPlotDict['title'] = ''
    myPlotDict['ylabel'] = r'Error in $H^{1}(\mathcal{D})$-norm [log]'
    myPlotDict['out_filename'] = dir_fig+'convg_lShaped_singular_interpolant.pdf'
    myPlotDict['xlim'] = [8E-3, 5E-1]
    # myPlotDict['ylim'] = [1E-4, 5]
    myPlotDict['ref_data_labels'] = ['$O(h^{%2.2f})$'%slope1, '$O(h^{%2.2f})$'%slope2]
    common.plotLogLogData([hMax, hMax], [err1, err2], [ref1, ref2], myPlotDict)


if __name__ == "__main__":

    show_plot = True
    save_plot = True
    
    dir_sol = '../../output/lShaped_interpolant/'
    dir_fig = '../../figures/lShaped_interpolant/'
   
    convg(dir_sol, dir_fig, show_plot, save_plot)
    
    
# End of file
