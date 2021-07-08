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


def convg(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot):

    if ptype == 1:
        fn = dir_sol+'convergenceFG_lShaped_test1_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'
    elif ptype == 2:
        fn = dir_sol+'convergenceFG_lShaped_test1_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
    hMax, ndofs, err1, err2 = load_data(fn)
    
    # linear fitting
    m = hMax.size
    linfit1 = np.polyfit(np.log(hMax[m-3:m]), np.log(err1[m-3:m]), 1)
    linfit2 = np.polyfit(np.log(hMax[m-3:m]), np.log(err2[m-3:m]), 1)
    ref1 = np.exp(np.polyval(linfit1, np.log(hMax))-0.5)
    ref2 = np.exp(np.polyval(linfit2, np.log(hMax))+0.5)
    slope1 = linfit1[0]
    slope2 = linfit2[0]
    
    myPlotDict = {}
    myPlotDict['show_plot'] = show_plot
    myPlotDict['save_plot'] = save_plot
    
    # plot error vs mesh size
    myPlotDict['xlabel'] = r'Meshsize $h$ [log]'
    myPlotDict['legend_loc'] = 'upper left'
    myPlotDict['data_markers'] = ['rs-', 'bo-']
    myPlotDict['data_labels'] = [r'$v$', r'$\mathbf{\sigma}$']
    myPlotDict['ylim'] = []
    myPlotDict['ref_data_markers'] = ['r--', 'b-.']
    
    if ptype == 1:
        myPlotDict['title'] = r'$p^{v}_{x} = p^{v}_{t} = p^{\mathbf{\sigma}}_{t} = p^{\mathbf{\sigma}}_{x} =$ '+str(deg)
    elif ptype == 2:
        myPlotDict['title'] = r'$p^{v}_{x} = p^{v}_{t} = p^{\mathbf{\sigma}}_{t} =$ '+str(deg)+r',   $p^{\mathbf{\sigma}}_{x} =$ '+str(deg-1)
    myPlotDict['ylabel'] = r'Error [log]'
    if ptype == 1:
        myPlotDict['out_filename'] = dir_fig+'convg_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
    elif ptype == 2:
        myPlotDict['out_filename'] = dir_fig+'convg_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.pdf'    
    myPlotDict['xlim'] = [8E-3, 5E-1]
    # myPlotDict['xlim'] = [5E-3, 5E-1]
    # myPlotDict['ylim'] = [1E-4, 5]
    myPlotDict['ref_data_labels'] = ['$O(h^{%2.2f})$'%slope1, '$O(h^{%2.2f})$'%slope2]
    common.plotLogLogData([hMax, hMax], [err1, err2], [ref1, ref2], myPlotDict)


def convgL2VsDG(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot):

    if ptype == 1:
        fn1 = dir_sol+'convergenceFG_lShaped_test1_L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceFG_lShaped_test1_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'    
    elif ptype == 2:
        fn1 = dir_sol+'convergenceFG_lShaped_test1_L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceFG_lShaped_test1_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
    hMax, ndofs, L2err1, L2err2 = load_data(fn1)
    hMax, ndofs, DGerr1, DGerr2 = load_data(fn2)
    
    errL2 = 0.5*(L2err1 + L2err2)
    errDg = np.sqrt(DGerr1**2 + DGerr2**2)
    print('Error L2: ', errL2)
    print('Error DG: ', errDg)
    
    # linear fitting
    m = hMax.size
    linfitL2 = np.polyfit(np.log(hMax[m-3:m]), np.log(errL2[m-3:m]), 1)
    linfitDg = np.polyfit(np.log(hMax[m-3:m]), np.log(errDg[m-3:m]), 1)
    refL2 = np.exp(np.polyval(linfitL2, np.log(hMax))-0.5)
    refDg = np.exp(np.polyval(linfitDg, np.log(hMax))+0.5)
    slopeL2 = linfitL2[0]
    slopeDg = linfitDg[0]
    
    myPlotDict = {}
    myPlotDict['show_plot'] = show_plot
    myPlotDict['save_plot'] = save_plot
    
    # plot error vs mesh size
    myPlotDict['xlabel'] = r'Meshsize $h$ [log]'
    myPlotDict['legend_loc'] = 'upper left'
    myPlotDict['data_markers'] = ['rs-', 'bo-']
    myPlotDict['data_labels'] = [r'$L^{2}(\mathcal{D}\times\{T\})$', 'DG(Q)']
    myPlotDict['ylim'] = []
    myPlotDict['ref_data_markers'] = ['r--', 'b-.']
    
    if ptype == 1:
        myPlotDict['title'] = r'$p^{v}_{x} = p^{v}_{t} = p^{\mathbf{\sigma}}_{t} = p^{\mathbf{\sigma}}_{x} =$ '+str(deg)
    elif ptype == 2:
        myPlotDict['title'] = r'$p^{v}_{x} = p^{v}_{t} = p^{\mathbf{\sigma}}_{t} =$ '+str(deg)+r',   $p^{\mathbf{\sigma}}_{x} =$ '+str(deg-1)
    myPlotDict['ylabel'] = r'Error [log]'
    if ptype == 1:
        myPlotDict['out_filename'] = dir_fig+'convg_L2vsDG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
    elif ptype == 2:
        myPlotDict['out_filename'] = dir_fig+'convg_L2vsDG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.pdf'    
    myPlotDict['xlim'] = [8E-3, 5E-1]
    # myPlotDict['ylim'] = [1E-4, 5]
    myPlotDict['ref_data_labels'] = ['$O(h^{%2.2f})$'%slopeL2, '$O(h^{%2.2f})$'%slopeDg]
    common.plotLogLogData([hMax, hMax], [errL2, errDg], [refL2, refDg], myPlotDict)


def convgFGvsSGWrtNdofsRelL2(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, plot_id):

    if ptype == 1:
        fn1 = dir_sol+'convergenceFG_lShaped_test1_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_lShaped_test1_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'    
    elif ptype == 2:
        fn1 = dir_sol+'convergenceFG_lShaped_test1_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_lShaped_test1_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
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
    myPlotDict['ylabel'] = r'Rel. $L^2(\mathcal{D}_{\mathbf{x}}\times\{T\})$ error [log]'
    if ptype == 1:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSG_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
    elif ptype == 2:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSG_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.pdf'    
    myPlotDict['xlim'] = [5E+3, 5E+8]
    # myPlotDict['ylim'] = [1E-4, 5]
    myPlotDict['ref_data_labels'] = ['$O({M_{L}}^{%2.2f})$'%slopeFG, '$O({M_{L}}^{%2.2f})$'%slopeSG]
    common.plotLogLogData([ndofsFG, ndofsSG], [errFG, errSG], [refFG, refSG], myPlotDict)


def convgFGvsSGWrtNdofsDG(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, plot_id):

    if ptype == 1:
        fn1 = dir_sol+'convergenceFG_lShaped_test1_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_lShaped_test1_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'    
    elif ptype == 2:
        fn1 = dir_sol+'convergenceFG_lShaped_test1_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_lShaped_test1_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
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
    myPlotDict['xlim'] = [5E+3, 5E+8]
    # myPlotDict['ylim'] = [1E-4, 5]
    myPlotDict['ref_data_labels'] = ['$O({M_{L}}^{%2.2f})$'%slopeFG, '$O({M_{L}}^{%2.2f})$'%slopeSG]
    common.plotLogLogData([ndofsFG, ndofsSG], [errFG, errSG], [refFG, refSG], myPlotDict)


if __name__ == "__main__":

    show_plot = True
    save_plot = True
    """
    dir_sol = '../../output/lShaped_test1/qu/'
    dir_fig = '../../figures/lShaped_test1/qu/'
   
    ptype = 1
    deg = 3
    stab = 1
    convg(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot)
    """
    """
    dir_sol = '../../output/lShaped_test1/br/'
    dir_fig = '../../figures/lShaped_test1/br/'
    
    ptype = 2
    deg = 3
    stab = 4
    convgL2VsDG(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot)
    """
    """
    dir_sol = '../../output/lShaped_test1/rg_nc/'
    dir_fig = '../../figures/lShaped_test1/rg_nc/'
   
    ptype = 1
    deg = 2
    stab = 1
    # convg(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot)
    convgFGvsSGWrtNdofsRelL2(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, 1)
    convgFGvsSGWrtNdofsDG(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, 1)
    """
    
    dir_sol = '../../output/lShaped_test1/br/'
    dir_fig = '../../figures/lShaped_test1/br/'
    
    ptype = 1
    deg = 2
    stab = 1
    convgFGvsSGWrtNdofsRelL2(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, 1)
    convgFGvsSGWrtNdofsDG(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, 1)
    
    
# End of file
