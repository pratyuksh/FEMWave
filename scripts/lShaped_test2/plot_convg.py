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

# L2L2 error
def convgFGvsSGWrtNdofsRelL2L2(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, plot_id):

    if ptype == 1:
        fn1 = dir_sol+'convergenceFG_lShaped_test2_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_lShaped_test2_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'    
    elif ptype == 2:
        fn1 = dir_sol+'convergenceFG_lShaped_test2_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_lShaped_test2_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
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
    myPlotDict['ylabel'] = r'Rel. $L^{2}(Q)$ error [log]'
    if ptype == 1:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSG_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
    elif ptype == 2:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSg_relL2L2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.pdf'    
    myPlotDict['xlim'] = [5E+3, 5E+8]
    # myPlotDict['ylim'] = [1E-4, 5]
    myPlotDict['ref_data_labels'] = ['$O({M_{L}}^{%2.2f})$'%slopeFG, '$O({M_{L}}^{%2.2f})$'%slopeSG]
    common.plotLogLogData([ndofsFG, ndofsSG], [errFG, errSG], [refFG, refSG], myPlotDict)

# DG error
def convgFGvsSGWrtNdofsDG(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, plot_id):

    if ptype == 1:
        fn1 = dir_sol+'convergenceFG_lShaped_test2_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_lShaped_test2_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'    
    elif ptype == 2:
        fn1 = dir_sol+'convergenceFG_lShaped_test2_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
        fn2 = dir_sol+'convergenceSG_lShaped_test2_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
    hMaxFG, ndofsFG, errFG1, errFG2 = load_data(fn1)
    hMaxSG, ndofsSG, errSG1, errSG2 = load_data(fn2)
    
    errFG = np.sqrt(errFG1**2 + errFG2**2)
    errSG = np.sqrt(errSG1**2 + errSG2**2)
    #errFG = errFG2
    #errSG = errSG2
    print('Error FG: ', errFG)
    print('Error SG: ', errSG)
    
    # linear fitting
    m1 = ndofsFG.size
    m2 = ndofsSG.size
    linfitFG = np.polyfit(np.log(ndofsFG[m1-3:m1]), np.log(errFG[m1-3:m1]), 1)
    linfitSG = np.polyfit(np.log(ndofsSG[m2-3:m2]), np.log(errSG[m2-3:m2]), 1)
    #linfitFG = np.polyfit(np.log(ndofsFG[m1-4:m1-1]), np.log(errFG[m1-4:m1-1]), 1)
    #linfitSG = np.polyfit(np.log(ndofsSG[m2-4:m2-1]), np.log(errSG[m2-4:m2-1]), 1)
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
    #myPlotDict['ylabel'] = r'DG error, space-like [log]'
    #myPlotDict['ylabel'] = r'DG error, time-like [log]'
    if ptype == 1:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSG_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
        #myPlotDict['out_filename'] = dir_fig+'convgFGvsSG_DGspace_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
        #myPlotDict['out_filename'] = dir_fig+'convgFGvsSG_DGtime_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
    elif ptype == 2:
        myPlotDict['out_filename'] = dir_fig+'convgFGvsSg_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.pdf'    
    # myPlotDict['xlim'] = [5E+3, 5E+8]
    myPlotDict['xlim'] = [1E+3, 5E+8]
    # myPlotDict['ylim'] = [1E-4, 5]
    myPlotDict['ref_data_labels'] = ['$O({M_{L}}^{%2.2f})$'%slopeFG, '$O({M_{L}}^{%2.2f})$'%slopeSG]
    common.plotLogLogData([ndofsFG, ndofsSG], [errFG, errSG], [refFG, refSG], myPlotDict)


# DG error, SG, components separately
def convgSGWrtNdofsDG(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, plot_id):

    if ptype == 1:
        fn1 = dir_sol+'convergenceSG_lShaped_test2_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'_p.json'
        fn2 = dir_sol+'convergenceSG_lShaped_test2_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'_v.json'    
    elif ptype == 2:
        fn1 = dir_sol+'convergenceSG_lShaped_test2_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'_p.json'
        fn2 = dir_sol+'convergenceSG_lShaped_test2_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'_v.json'
    hMaxSGp, ndofsSGp, errSGp1, errSGp2 = load_data(fn1)
    hMaxSGv, ndofsSGv, errSGv1, errSGv2 = load_data(fn2)
    
    print('Error p: ', errSGp1, ' ', errSGp2)
    print('Error v: ', errSGv1, ' ', errSGv2)
    
    if plot_id == 1:
    
        # linear fitting
        m1 = ndofsSGp.size
        m2 = ndofsSGv.size
        linfitSGp = np.polyfit(np.log(ndofsSGp[m1-3:m1]), np.log(errSGp1[m1-3:m1]), 1)
        linfitSGv = np.polyfit(np.log(ndofsSGv[m2-3:m2]), np.log(errSGv1[m2-3:m2]), 1)
        refSGp = np.exp(np.polyval(linfitSGp, np.log(ndofsSGp))+0.5)
        refSGv = np.exp(np.polyval(linfitSGv, np.log(ndofsSGv))-0.5)
        slopeSGp = linfitSGp[0]
        slopeSGv = linfitSGv[0]
    
        myPlotDict = {}
        myPlotDict['show_plot'] = show_plot
        myPlotDict['save_plot'] = save_plot
    
        # plot error vs mesh size
        myPlotDict['xlabel'] = r'Number of degrees of freedom, $M_{L}$ [log]'
        myPlotDict['legend_loc'] = 'upper right'
        myPlotDict['data_markers'] = ['rs-', 'bo-']
        myPlotDict['data_labels'] = [r'p='+str(deg)+', SG, $v$', r'p='+str(deg)+', SG, $\mathbf{\sigma}$']
        myPlotDict['ylim'] = []
        myPlotDict['ref_data_markers'] = ['r--', 'b-.']
    
        myPlotDict['title'] = ''
        myPlotDict['ylabel'] = r'DG error, space-like [log]'
        if ptype == 1:
            myPlotDict['out_filename'] = dir_fig+'convgSG_DGspace_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
        elif ptype == 2:
            myPlotDict['out_filename'] = dir_fig+'convgSG_DGspace_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.pdf'    
        myPlotDict['xlim'] = [1E+5, 5E+8]
        # myPlotDict['ylim'] = [1E-4, 5]
        myPlotDict['ref_data_labels'] = ['$O({M_{L}}^{%2.2f})$'%slopeSGp, '$O({M_{L}}^{%2.2f})$'%slopeSGv]
        common.plotLogLogData([ndofsSGp, ndofsSGv], [errSGp1, errSGv1], [refSGp, refSGv], myPlotDict)
    
    elif plot_id == 2:
    
        # linear fitting
        m1 = ndofsSGp.size
        m2 = ndofsSGv.size
        #linfitSGp = np.polyfit(np.log(ndofsSGp[m1-3:m1]), np.log(errSGp2[m1-3:m1]), 1)
        #linfitSGv = np.polyfit(np.log(ndofsSGv[m2-3:m2]), np.log(errSGv2[m2-3:m2]), 1)
        linfitSGp = np.polyfit(np.log(ndofsSGp[m1-3:m1-1]), np.log(errSGp2[m1-3:m1-1]), 1)
        linfitSGv = np.polyfit(np.log(ndofsSGv[m2-3:m2-1]), np.log(errSGv2[m2-3:m2-1]), 1)
        refSGp = np.exp(np.polyval(linfitSGp, np.log(ndofsSGp))-0.5)
        refSGv = np.exp(np.polyval(linfitSGv, np.log(ndofsSGv))+0.5)
        slopeSGp = linfitSGp[0]
        slopeSGv = linfitSGv[0]
    
        myPlotDict = {}
        myPlotDict['show_plot'] = show_plot
        myPlotDict['save_plot'] = save_plot
    
        # plot error vs mesh size
        myPlotDict['xlabel'] = r'Number of degrees of freedom, $M_{L}$ [log]'
        myPlotDict['legend_loc'] = 'upper right'
        myPlotDict['data_markers'] = ['rs-', 'bo-']
        myPlotDict['data_labels'] = [r'p='+str(deg)+', SG, $v$', r'p='+str(deg)+', SG, $\mathbf{\sigma}$']
        myPlotDict['ylim'] = []
        myPlotDict['ref_data_markers'] = ['r--', 'b-.']
    
        myPlotDict['title'] = ''
        myPlotDict['ylabel'] = r'DG error, time-like [log]'
        if ptype == 1:
            myPlotDict['out_filename'] = dir_fig+'convgSG_DGtime_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.pdf'
        elif ptype == 2:
            myPlotDict['out_filename'] = dir_fig+'convgSG_DGtime_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.pdf'    
        myPlotDict['xlim'] = [1E+5, 5E+8]
        # myPlotDict['ylim'] = [1E-4, 5]
        myPlotDict['ref_data_labels'] = ['$O({M_{L}}^{%2.2f})$'%slopeSGp, '$O({M_{L}}^{%2.2f})$'%slopeSGv]
        common.plotLogLogData([ndofsSGp, ndofsSGv], [errSGp2, errSGv2], [refSGp, refSGv], myPlotDict)

if __name__ == "__main__":

    show_plot = True
    save_plot = True
    
    #dir_sol = '../../output/lShaped_test2/br/projection/'
    #dir_fig = '../../figures/lShaped_test2/br/projection/'
    
    #dir_sol = '../../output/lShaped_test2/br/'
    #dir_fig = '../../figures/lShaped_test2/br/'
    
    dir_sol = '../../output/lShaped_test2/rg_nc/'
    dir_fig = '../../figures/lShaped_test2/rg_nc/'
    
    ptype = 1
    deg = 2
    stab = 1
    convgFGvsSGWrtNdofsRelL2L2(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, 1)
    # convgFGvsSGWrtNdofsDG(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, 1)
    # convgSGWrtNdofsDG(ptype, stab, deg, dir_sol, dir_fig, show_plot, save_plot, 2)
    
# End of file
