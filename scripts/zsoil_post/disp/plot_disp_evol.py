# -*- coding: cp1252 -*-
# @description Plotting evolution of nodal displacements.
# @output displacement plot
# @author Matthias Preisig
# @date 2017/08/18
import os,math
import numpy as np
import matplotlib.pyplot as plt
import pickle

from zsoil_tools import zsoil_results as zr

pathname = '..'
pblist = ['M1352_essai_stat_CPTUKF2_CBMCMW_v1_1n']

for kp,prob in enumerate(pblist):

    figname = prob
    try:
        res = pickle.load(open(prob+'_disp.p', "rb" ))
        print(prob+' loaded')
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

    crds = [np.array([0,0.8,0])]
    volnds = set()
    for ke in range(res.nVolumics):
        if res.vol.mat[ke]==11:
            inel = res.vol.inel[ke]
            for kn in inel:
                volnds.add(kn-1)
    volnds = list(volnds)
    nlist = [-1 for k in crds]
    for kn in range(res.nNodes):
        if kn in volnds:
            for kk,crd in enumerate(crds):
                if abs(res.coords[0][kn]-crd[0])<1e-3:
                    if abs(res.coords[1][kn]-crd[1])<1e-3:
                        if abs(res.coords[2][kn]-crd[2])<1e-3:
                            nlist[kk] = kn
                            break
    step0 = res.steps[0]

    tx = []
    t1 = 300
    t2 = 1200
    uv = [[] for kn in nlist]
    for kt in res.out_steps:
        step = res.steps[kt]
        if step.time>t2:
            tx.append(step.time-t2+t1+100)
        elif step.time>t1:
            tx.append(t1+(step.time-t1)/9)
        else:
            tx.append(step.time)
        for kk,kn in enumerate(nlist):
            uv[kk].append(1e3*(step.nodal.disp[1][kn]-
                               step0.nodal.disp[1][kn]))
    
    txticks0 = [120,180,220,240,260,280,300,1270,1280,1300,
                1320,1340,1360,1380,1400,1420,1440,1460]
    xticklabels = ['colonnes +\ninfrastructure',
                   '120 kN (poids mort)',''
                   '200 kN','300 kN','400 kN','500 kN','120 kN','',
                   '200 kN','300 kN','400 kN','500 kN','600 kN', 
                   '700 kN','800 kN','900 kN','1000 kN','']                 
    txticks = []
    for t in txticks0:
        if t>t2:
            txticks.append(t-t2+t1+100)
        elif t>t1:
            txticks.append(t1+(t-t1)/9)
        else:
            txticks.append(t)

    lstr = ['bs-','r*-','g<-','mo-']
    lstr = ['b','g','r','m','y','c','b--','g--','r--','m--','y--','c--']

    fig = plt.figure(figsize=(12,8))
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    fig.text(0.01,0.01,prob,size=8)
        
    ax = fig.add_subplot(111)

    for kk,u in enumerate(uv):
        ax.plot(tx,u,lstr[kk])

    ax.set_ylabel(u'Tassement [mm]',size=14)
    ax.legend()
    ax.set_xticks(txticks)
    ax.set_xticklabels(xticklabels)
    plt.setp(ax.xaxis.get_majorticklabels(),rotation=30,ha='right')
    ax.set_xlim([0,600])
    ax.grid('on')
                
    fig.tight_layout()
    fig.savefig(prob+'_dv')
    plt.close(fig)

        

            




