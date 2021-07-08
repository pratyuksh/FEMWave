#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import numpy as np
from decimal import Decimal


def load_json(filename):
    with open(filename, "r") as f:
        return json.load(f)

def load_data(filename):
    j = load_json(filename)
    return np.array(j["h_max"]), np.array(j["ndofs"], dtype=np.int32), np.array(j["pressure"]), np.array(j["velocity"])


def convgRateWrtH (data_x, data_y):
    N = data_x.size
    #z = np.polyfit(np.log2(data_x), np.log2(data_y), 1)
    #z = np.polyfit(np.log2(data_x[1:N]), np.log2(data_y[1:N]), 1)
    z = np.polyfit(np.log2(data_x[N-3:N]), np.log2(data_y[N-3:N]), 1)
    return round(z[0],2)

# Relative L2-error
def evalRateErrRelL2FG(dir_sol, ptype, stab, deg):
    
    if ptype == 1:
        fn = dir_sol+'convergence_unit_square_test1_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'    
    elif ptype == 2:
        fn = dir_sol+'convergence_unit_square_test1_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
    hMax, ndofs, errp, errv = load_data(fn)
    
    print("\n\n")
    for k in range(errp.shape[0]-1,errp.shape[0]):
        print(k, "%.4E"%Decimal(errp[k]), "%.4E"%Decimal(errv[k]))
    print("\n")
    for k in range(errp.size-3,errp.size-1):
        val1 = np.log2(errp[k]) - np.log2(errp[k+1])
        val2 = np.log2(errv[k]) - np.log2(errv[k+1])
        print(k, "%4.2f"%val1, "%4.2f"%val2)
    
    ratep = convgRateWrtH (hMax, errp[:])
    ratev = convgRateWrtH (hMax, errv[:])
    print('\nConvergence rate of pressure: ', ratep)
    print('Convergence rate of velocity: ', ratev)

# L2-error
def evalRateErrL2FG(dir_sol, ptype, stab, deg):
    
    if ptype == 1:
        fn = dir_sol+'convergence_unit_square_test1_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'    
    elif ptype == 2:
        fn = dir_sol+'convergence_unit_square_test1_relL2_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
    hMax, ndofs, err1, err2 = load_data(fn)
    
    errL2 = 0.5*(err1 + err2)
    
    print("\n\n")
    for k in range(errL2.shape[0]-1, errL2.shape[0]):
        print(k, "%.4E"%Decimal(errL2[k]))
    print("\n")
    for k in range(errL2.size-3, errL2.size-1):
        val = np.log2(errL2[k]) - np.log2(errL2[k+1])
        print(k, "%4.2f"%val)
    
    rate = convgRateWrtH (hMax, errL2[:])
    print('\nConvergence rate of L2 error: ', rate)

# DG-error
def evalRateErrDgFG(dir_sol, ptype, stab, deg):
    
    if ptype == 1:
        fn = dir_sol+'convergence_unit_square_test1_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg)+'_stab'+str(stab)+'.json'    
    elif ptype == 2:
        fn = dir_sol+'convergence_unit_square_test1_DG_degt'+str(deg)+'_deg1x'+str(deg)+'_deg2x'+str(deg-1)+'_stab'+str(stab)+'.json'
    hMax, ndofs, err1, err2 = load_data(fn)
    
    errDg = np.sqrt(err1**2 + err2**2)
    
    print("\n\n")
    for k in range(errDg.shape[0]-1, errDg.shape[0]):
        print(k, "%.4E"%Decimal(errDg[k]))
    print("\n")
    for k in range(errDg.size-3,errDg.size-1):
        val = np.log2(errDg[k]) - np.log2(errDg[k+1])
        print(k, "%4.2f"%val)
    
    rate = convgRateWrtH (hMax, errDg[:])
    print('\nConvergence rate of Dg error: ', rate)


if __name__=='__main__':
    
    dir_sol = '../../output/unitSquare_test1/qu/'
    
    print("\n\nFull-grid smooth solution:\n")
    ptype = 1
    deg = 1
    for stabParams in range(1,5):
        #evalRateErrRelL2FG(dir_sol, ptype, stabParams, deg)
        evalRateErrL2FG(dir_sol, ptype, stabParams, deg)
        #evalRateErrDgFG(dir_sol, ptype, stabParams, deg)
    

# End of file
