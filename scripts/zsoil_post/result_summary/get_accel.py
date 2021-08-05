# -*- coding: cp1252 -*-
# @description Extracting time series from zsoil results
# @input zsoil results
# @output time history in txt-file
# @author Matthias Preisig
# @date 2016/05/31
import numpy,os,math
from numpy import linalg as la
from numpy import fft
import matplotlib.pyplot as plt
from scipy import interpolate


import sys
sys.path.append('D:/Mandats/python/zsoil_tools')
sys.path.append('D:/Mandats/python/utils/response_spectrum')
import zsoil_results as zr
import resp_spect
import time

def ag(t):
    if t>max(tx):
        return 0
    else:
        return finterp(t)

def get_beam_node(res,x,y,mat):
    for ke in range(res.nBeams):
        if res.beam.mat[ke]==mat:
            for kn in res.beam.inel[ke]:
                if abs(res.coords[0][kn-1]-x)<1e-6 and abs(res.coords[1][kn-1]-y)<1e-6:
                    return kn-1
    return -1

def get_vol_node(res,x,y):
    nlist = []
    for kn in range(res.nNodes):
        if abs(res.coords[0][kn]-x)<1e-6 and abs(res.coords[1][kn]-y)<1e-6:
            nlist.append(kn+1)        
    for ke in range(res.nVolumics):
        for kn in nlist:
            if kn in res.vol.inel[ke]:
                xm = 0.25*sum([res.coords[0][kn-1] for kn in res.vol.inel[ke]])
                ym = 0.25*sum([res.coords[1][kn-1] for kn in res.vol.inel[ke]])
                if xm<x and ym<y:
                    return kn-1
    return -1

pathname = '//CALCUL-APR11/Mandats Disc C/M890_MurGareBern/2015'
    

figname = 'cdb'
pblist = ['schnitt1_2Dequiv_Projekt_acc01',
          'schnitt1_2Dequiv_fBase_Projekt_acc01',
          'schnitt1_2Dequiv_fBaseRA_Projekt_acc01']
labstr = [u'Dämpfer','Fixe Basis','Fix, Rayleigh']
ndlist = [1988,1989,1990]
nd0 = 1200

figname = 'x1p3_gf14_ohnePGS'
pblist = ['schnitt1_2Dequiv_Projekt_ohnePGS_acc01x1p3_gf14',
          'schnitt1_2Dequiv_Projekt_ohnePGS_acc142yx1p3_gf14',
          'schnitt1_2Dequiv_Projekt_ohnePGS_acc369xx12x1p3_gf14']
labstr = ['acc01','acc142y','acc369xx12']

figname = 'Schnitt1_x1p3_gf14'
pblist = ['schnitt1_2Dequiv_Projekt_acc01x1p3_gf14',
          'schnitt1_2Dequiv_Projekt_acc142yx1p3_gf14',
          'schnitt1_2Dequiv_Projekt_acc369xx12x1p3_gf14']
ndlist = [173]# ZS numbers

##figname = 'Schnitt2_x1p3_gf14_05xP0'
##pblist = ['schnitt2_2Dequiv_Projekt_acc01x1p3_gf14_05xP0',
##          'schnitt2_2Dequiv_Projekt_acc142yx1p3_gf14_05xP0',
##          'schnitt2_2Dequiv_Projekt_acc369xx12x1p3_gf14_05xP0',
##          'schnitt2_2Dequiv_Projekt_ohnePGS_acc01x1p3_gf14_05xP0',
##          'schnitt2_2Dequiv_Projekt_ohnePGS_acc142yx1p3_gf14_05xP0',
##          'schnitt2_2Dequiv_Projekt_ohnePGS_acc369xx12x1p3_gf14_05xP0']
##spacing = 2.95
##ktruss = [1,2,3]
##blist = [19,20,21]
##bbase = 30
##ndlist = [1419,1420,1421]
##nd0 = 802

for kf,prob in enumerate(pblist):
    res = zr.zsoil_results(pathname,prob)
    res.read_rcf()
    res.read_his()
    tsteps = range(0,len(res.steps))
    tsteps = []
    for kt,step in enumerate(res.steps):
        if step.conv_status==-1 and step.time>=3:
            tsteps.append(kt)
    res.out_steps = tsteps
    res.read_dat()
    res.read_LTF()
    res.read_s00('/v500',res_type='accelerations')

    finterp = interpolate.interp1d(res.LTF[4][0],res.LTF[4][1])

    tx = []
    dispx = [[] for k in ndlist]
    dx0 = []
    for kkt,kt in enumerate(res.out_steps):
        step = res.steps[kt]
        t = step.time
        tx.append(t-3)
        xbase = ag(t)*1.4
##        xbase = res.steps[kt].nodal.disp[0][nd0-1]
        for kkn,kn in enumerate(ndlist):
            dx = step.nodal.a_disp[0][kn-1]
            dispx[kkn].append(dx+xbase)

    
    for kl in range(len(dispx)):
        of = open(prob+'_acc_terrain_node%i.txt'%(ndlist[kl]),'w')
        for kt in range(len(tx)):
            of.write('%1.4e %1.4e\n'%(tx[kt],dispx[kl][kt]))
        of.close()





