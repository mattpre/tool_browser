# -*- coding: cp1252 -*-
# @description Plotting evolution of nodal displacements.
# @output displacement plot
# @author Matthias Preisig
# @date 2017/08/18
import os,math
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle

from zsoil_tools import zsoil_results as zr

pathname = '..'
pblist = ['M1133_full_v1_4_l=100%','M1133_full_v2_l=100%']
pblist = ['M1133_Full_v2_l=100%_E10_VoussMet']

for kp,prob in enumerate(pblist):

    figname = prob
    try:
        res = pickle.load(open(prob+'_disp.p', "rb" ))
        print prob+' loaded'
        tsteps = res.out_steps
        tvect = [res.steps[kt].time for kt in tsteps]
    except:
        res = zr(pathname,prob)
        res.read_rcf()
        res.read_his()

        tsteps = []
        tsteps_plot = []
        for kt,step in enumerate(res.steps):
            if step.sf==0 and step.conv_status==-1:
                tsteps.append(kt)

        res.out_steps = tsteps
        res.read_dat()
        res.read_s00()
        pickle.dump(res, open(prob+'_disp.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    hcrds = [np.array([1.68010e+01,2.23000e+01,-2.84842e+01]),
             np.array([3.13969e+01,2.23000e+01,-1.03852e+01]),
             np.array([3.19070e+01,2.09250e+01,8.69259e+00]),
             np.array([0.00000e+00,2.23000e+01,3.30700e+01]),
             np.array([-3.19070e+01,2.09250e+01,8.69259e+00]),
             np.array([-3.13969e+01,2.09250e+01,-1.03852e+01]),
             np.array([0.00000e+00,2.23000e+01,-3.30700e+01])]
    hlabs = ['Nord de L17N, Est',
             'Entre L16 et L17N, Est',
             'Entre L16 et L17S, Est',
             'Sud de L17S, sur l\'axe',
             'Entre L16 et L17S, Ouest',
             'Entre L16 et L17N, Ouest',
             'Nord de L17N, sur l\'axe']
    vcrds = [np.array([4.04606e+01,5.00000e+01,-2.17519e+01]),
             np.array([2.06696e+00,4.40000e+01,-3.68605e+01]),
             np.array([-3.45153e+01,4.40000e+01,7.53607e-15]),
             np.array([0.00000e+00,4.40000e+01,3.73025e+01]),
             np.array([3.45153e+01,4.40000e+01,7.53607e-15]),
             np.array([3.62460e+01,4.20000e+01,-2.23071e+01]),
             np.array([4.83033e+01,4.20000e+01,-3.77146e+00]),
             np.array([6.88199e+01,4.20000e+01,2.92640e+01]),
             np.array([8.05865e+01,4.20000e+01,4.74675e+01])]
    vlabs = ['Digue GC',
             'Nord',
             'Ouest',
             'Sud',
             'Est',
             u'GC, culée Nord',
             u'GC, appui Nord',
             u'GC, appui Sud',
             u'GC, culée Sud']
    shnds = set()
    for ke in range(res.nShells):
        if res.shell.mat[ke]==11:
            inel = res.shell.inel[ke]
            for kn in inel:
                shnds.add(kn-1)
    shnds = list(shnds)
    nhlist = [-1 for k in hlabs]
    nvlist = [-1 for k in vlabs]    
    for kn in range(res.nNodes):
        if kn in shnds:
            for kk,crd in enumerate(hcrds):
                if abs(res.coords[0][kn]-crd[0])<1e-3:
                    if abs(res.coords[1][kn]-crd[1])<1e-3:
                        if abs(res.coords[2][kn]-crd[2])<1e-3:
                            nhlist[kk] = kn
                            break
        for kk,crd in enumerate(vcrds):
            if abs(res.coords[0][kn]-crd[0])<1e-3:
                if abs(res.coords[1][kn]-crd[1])<1e-3:
                    if abs(res.coords[2][kn]-crd[2])<1e-3:
                        nvlist[kk] = kn
                        break
    step0 = res.steps[tsteps[2]]

    tx = []
    uh = [[] for kn in nhlist]
    uv = [[] for kn in nvlist]
    for kt in tsteps[2:]:
        step = res.steps[kt]
        tx.append(step.time)
        for kk,kn in enumerate(nhlist):
            uh[kk].append(1e3*((step.nodal.disp[0][kn]**2+step.nodal.disp[2][kn]**2)**0.5-
                               (step0.nodal.disp[0][kn]**2+step0.nodal.disp[2][kn]**2)**0.5))
        for kk,kn in enumerate(nvlist):
            uv[kk].append(1e3*(step.nodal.disp[1][kn]-step0.nodal.disp[1][kn]))

    lstr = ['bs-','r*-','g<-','mo-']
    lstr = ['b','g','r','m','y','c','b--','g--','r--','m--','y--','c--']

    fig = plt.figure(figsize=(12,8))
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    fig.text(0.01,0.01,prob,size=8)
        
    ax = fig.add_subplot(111)

    for kk,u in enumerate(uh):
        ax.plot(tx,u,lstr[kk],label=hlabs[kk])

    ax.set_ylabel(u'Déplacement horizontal [mm]',size=14)
    ax.legend()
    ax.set_xticks([12,14,25,35,37,48,51,62,64])
    etapes = [#u'réalisation pm',
              u'percement Est L16',
              u'percement Ouest L16',
              u'excavation S1',
              u'percement Est L17N',
              u'percement Ouest L17N',
              u'excavation S2',
              u'excavation f.f.',
              u'percement Est L17S',
              u'percement Ouest L17S']
    ax.set_xticklabels(etapes)
    plt.setp(ax.xaxis.get_majorticklabels(),rotation=30,ha='right')
    ax.set_xlim([0,75])
    ax.grid('on')
                
    fig.tight_layout()
    fig.savefig(prob+'_dh')
    plt.close(fig)

    fig = plt.figure(figsize=(12,8))
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    fig.text(0.01,0.01,prob,size=8)
        
    ax = fig.add_subplot(111)

    for kk,u in enumerate(uv):
        ax.plot(tx,u,lstr[kk],label=vlabs[kk])

    ax.set_ylabel(u'Tassement [mm]',size=14)
    ax.legend()
    ax.set_xticks([12,14,25,35,37,48,51,62,64])
    ax.set_xticklabels(etapes)
    plt.setp(ax.xaxis.get_majorticklabels(),rotation=30,ha='right')
    ax.set_xlim([0,75])
    ax.grid('on')
                
    fig.tight_layout()
    fig.savefig(prob+'_dv')
    plt.close(fig)

        

            




