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
    return np.array(j["h_max"]), np.array(j["ndofs"], dtype=np.int32), np.array(j["pressure"]), np.array(j["velocity"])


def convgFGvsSGWrtNdofsRelL2L2(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, plot_id):

    if ptype == 1:
        fn1 = dir_sol+'convergenceFG_unit_square_test4_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_unit_square_test4_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'    
    elif ptype == 2:
        fn1 = dir_sol+'convergenceFG_unit_square_test4_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_unit_square_test4_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
    hMaxFG, ndofsFG, errFG1, errFG2 = load_data(fn1)
    hMaxSG, ndofsSG, errSG1, errSG2 = load_data(fn2)
    
    errFG = np.sqrt(errFG1**2 + errFG2**2)
    errSG = np.sqrt(errSG1**2 + errSG2**2)
    print('Error FG: ', errFG)
    print('Error SG: ', errSG)
    
    # linear fitting
    m1 = ndofsFG.size
    m2 = ndofsSG.size
    linfitFG = np.polyfit(np.log(ndofsFG[m1-3:m1]), np.log(errFG[m1-3:m1]), 1)
    linfitSG = np.polyfit(np.log(ndofsSG[m2-3:m2]), np.log(errSG[m2-3:m2]), 1)
    #linfitFG = np.polyfit(np.log(ndofsFG), np.log(errFG), 1)
    #linfitSG = np.polyfit(np.log(ndofsSG), np.log(errSG), 1)
    refFG = np.exp(np.polyval(linfitFG, np.log(ndofsFG))+0.5)
    refSG = np.exp(np.polyval(linfitSG, np.log(ndofsSG))-0.5)
    slopeFG = linfitFG[0]
    slopeSG = linfitSG[0]
    
    myPlotDict = {}
    myPlotDict['show_plot'] = show_plot
    myPlotDict['save_plot'] = save_plot
    
    # plot error vs mesh size
    myPlotDict['xlabel'] = r'Number of degrees of freedom, $M_{L}$ [log]'
    myPlotDict['legend_loc'] = 'upper right'
    myPlotDict['data_markers'] = ['rs-', 'bo-']
    myPlotDict['data_labels'] = [r'p='+str(deg)+', FG', r'p='+str(deg)+', SG']
    myPlotDict['ylim'] = []
    myPlotDict['ref_data_markers'] = ['r--', 'b-.']
    
    if ptype == 1:
        myPlotDict['title'] = r'$p^{v}_{x} = p^{v}_{t} = p^{\mathbf{\sigma}}_{t} = p^{\mathbf{\sigma}}_{x} =$ '+str(deg)
    elif ptype == 2:
        myPlotDict['title'] = r'$p^{v}_{x} = p^{v}_{t} = p^{\mathbf{\sigma}}_{t} =$ '+str(deg)+r',   $p^{\mathbf{\sigma}}_{x} =$ '+str(deg-1)
    myPlotDict['ylabel'] = r'Rel. $L^2(Q)$ error [log]'
    if ptype == 1:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSG_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
    elif ptype == 2:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSg_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.pdf'    
    myPlotDict['xlim'] = [8E+1, 1E+8]
    # myPlotDict['ylim'] = [1E-4, 5]
    myPlotDict['ref_data_labels'] = ['$O({M_{L}}^{%2.2f})$'%slopeFG, '$O({M_{L}}^{%2.2f})$'%slopeSG]
    common.plotLogLogData([ndofsFG, ndofsSG], [errFG, errSG], [refFG, refSG], myPlotDict)


def convgFGvsSGWrtNdofsDG(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, plot_id):

    if ptype == 1:
        fn1 = dir_sol+'convergenceFG_unit_square_test4_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_unit_square_test4_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'    
    elif ptype == 2:
        fn1 = dir_sol+'convergenceFG_unit_square_test4_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_unit_square_test4_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
    hMaxFG, ndofsFG, errFG1, errFG2 = load_data(fn1)
    hMaxSG, ndofsSG, errSG1, errSG2 = load_data(fn2)
    
    errFG = np.sqrt(errFG1**2 + errFG2**2)
    errSG = np.sqrt(errSG1**2 + errSG2**2)
    print('Error FG: ', errFG)
    print('Error SG: ', errSG)
    
    # linear fitting
    m1 = ndofsFG.size
    m2 = ndofsSG.size
    linfitFG = np.polyfit(np.log(ndofsFG[m1-3:m1]), np.log(errFG[m1-3:m1]), 1)
    linfitSG = np.polyfit(np.log(ndofsSG[m2-3:m2]), np.log(errSG[m2-3:m2]), 1)
    #linfitFG = np.polyfit(np.log(ndofsFG), np.log(errFG), 1)
    #linfitSG = np.polyfit(np.log(ndofsSG), np.log(errSG), 1)
    refFG = np.exp(np.polyval(linfitFG, np.log(ndofsFG))+0.5)
    refSG = np.exp(np.polyval(linfitSG, np.log(ndofsSG))-0.5)
    slopeFG = linfitFG[0]
    slopeSG = linfitSG[0]
    
    myPlotDict = {}
    myPlotDict['show_plot'] = show_plot
    myPlotDict['save_plot'] = save_plot
    
    # plot error vs mesh size
    myPlotDict['xlabel'] = r'Number of degrees of freedom, $M_{L}$ [log]'
    myPlotDict['legend_loc'] = 'upper right'
    myPlotDict['data_markers'] = ['rs-', 'bo-']
    myPlotDict['data_labels'] = [r'p='+str(deg)+', FG', r'p='+str(deg)+', SG']
    myPlotDict['ylim'] = []
    myPlotDict['ref_data_markers'] = ['r--', 'b-.']
    
    if ptype == 1:
        myPlotDict['title'] = r'$p^{v}_{x} = p^{v}_{t} = p^{\mathbf{\sigma}}_{t} = p^{\mathbf{\sigma}}_{x} =$ '+str(deg)
    elif ptype == 2:
        myPlotDict['title'] = r'$p^{v}_{x} = p^{v}_{t} = p^{\mathbf{\sigma}}_{t} =$ '+str(deg)+r',   $p^{\mathbf{\sigma}}_{x} =$ '+str(deg-1)
    myPlotDict['ylabel'] = r'DG error [log]'
    if ptype == 1:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSG_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
    elif ptype == 2:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSg_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.pdf'    
    myPlotDict['xlim'] = [8E+1, 1E+8]
    # myPlotDict['ylim'] = [1E-4, 5]
    myPlotDict['ref_data_labels'] = ['$O({M_{L}}^{%2.2f})$'%slopeFG, '$O({M_{L}}^{%2.2f})$'%slopeSG]
    common.plotLogLogData([ndofsFG, ndofsSG], [errFG, errSG], [refFG, refSG], myPlotDict)


def convgFGvsSGWrtNdofsDGPlus(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, plot_id):

    if ptype == 1:
        fn1 = dir_sol+'convergenceFG_unit_square_test4_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceFG_unit_square_test4_DGPlus_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'
        fn3 = dir_sol+'convergenceSG_unit_square_test4_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'
        fn4 = dir_sol+'convergenceSG_unit_square_test4_DGPlus_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'    
    elif ptype == 2:
        fn1 = dir_sol+'convergenceFG_unit_square_test4_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceFG_unit_square_test4_DGPlus_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
        fn3 = dir_sol+'convergenceSG_unit_square_test4_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
        fn4 = dir_sol+'convergenceSG_unit_square_test4_DGPlus_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
    hMaxFG, ndofsFG, errFG11, errFG12 = load_data(fn1)
    hMaxFG, ndofsFG, errFG21, errFG22 = load_data(fn2)
    hMaxSG, ndofsSG, errSG11, errSG12 = load_data(fn3)
    hMaxSG, ndofsSG, errSG21, errSG22 = load_data(fn4)
    
    errFG = np.sqrt(errFG11**2 + errFG12**2 + errFG21**2 + errFG22**2)
    errSG = np.sqrt(errSG11**2 + errSG12**2 + errSG21**2 + errSG22**2)
    print('Error FG: ', errFG)
    print('Error SG: ', errSG)
    
    # linear fitting
    m1 = ndofsFG.size
    m2 = ndofsSG.size
    linfitFG = np.polyfit(np.log(ndofsFG[m1-3:m1]), np.log(errFG[m1-3:m1]), 1)
    linfitSG = np.polyfit(np.log(ndofsSG[m2-3:m2]), np.log(errSG[m2-3:m2]), 1)
    #linfitFG = np.polyfit(np.log(ndofsFG), np.log(errFG), 1)
    #linfitSG = np.polyfit(np.log(ndofsSG), np.log(errSG), 1)
    refFG = np.exp(np.polyval(linfitFG, np.log(ndofsFG))+0.5)
    refSG = np.exp(np.polyval(linfitSG, np.log(ndofsSG))-0.5)
    slopeFG = linfitFG[0]
    slopeSG = linfitSG[0]
    
    myPlotDict = {}
    myPlotDict['show_plot'] = show_plot
    myPlotDict['save_plot'] = save_plot
    
    # plot error vs mesh size
    myPlotDict['xlabel'] = r'Number of degrees of freedom, $M_{L}$ [log]'
    myPlotDict['legend_loc'] = 'upper right'
    myPlotDict['data_markers'] = ['rs-', 'bo-']
    myPlotDict['data_labels'] = [r'p='+str(deg)+', FG', r'p='+str(deg)+', SG']
    myPlotDict['ylim'] = []
    myPlotDict['ref_data_markers'] = ['r--', 'b-.']
    
    if ptype == 1:
        myPlotDict['title'] = r'$p^{v}_{x} = p^{v}_{t} = p^{\mathbf{\sigma}}_{t} = p^{\mathbf{\sigma}}_{x} =$ '+str(deg)
    elif ptype == 2:
        myPlotDict['title'] = r'$p^{v}_{x} = p^{v}_{t} = p^{\mathbf{\sigma}}_{t} =$ '+str(deg)+r',   $p^{\mathbf{\sigma}}_{x} =$ '+str(deg-1)
    myPlotDict['ylabel'] = r'DG error [log]'
    if ptype == 1:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSG_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
    elif ptype == 2:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSg_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.pdf'    
    myPlotDict['xlim'] = [8E+1, 1E+8]
    # myPlotDict['ylim'] = [1E-4, 5]
    myPlotDict['ref_data_labels'] = ['$O({M_{L}}^{%2.2f})$'%slopeFG, '$O({M_{L}}^{%2.2f})$'%slopeSG]
    common.plotLogLogData([ndofsFG, ndofsSG], [errFG, errSG], [refFG, refSG], myPlotDict)


if __name__ == "__main__":

    show_plot = True
    save_plot = True
    
    dir_sol = '../../output/unitSquare_test4/qu/'
    dir_fig = '../../figures/unitSquare_test4/'
    
    #dir_sol = '../../output/unitSquare_test4/qu/projection/'
    #dir_fig = '../../figures/unitSquare_test4/projection/'
    
    ptype = 1
    deg = 2
    stab = 1
    convgFGvsSGWrtNdofsRelL2L2(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, 1)
    # convgFGvsSGWrtNdofsDG(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, 1)
    # convgFGvsSGWrtNdofsDGPlus(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, 1)
    
    
# End of file
