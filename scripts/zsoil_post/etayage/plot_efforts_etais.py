# -*- coding: cp1252 -*-
# @description Plotting struts with normal forces or section data.
# @input zsoil results
# @output png 
# @author Matthias Preisig
# @date 2018/31/05
import os,math,cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import cPickle as pickle


from zsoil_tools import zsoil_results as zr


pathname = '..'
pblist = ['M1119_3D_v4','M1119_3D_v4_ROR']
pblist = ['M1119_3D_v4_ROR_P02']


for kf,prob in enumerate(pblist):
    figname = prob
    try:
        res = pickle.load(open(prob+'_etais.p', "rb" ))
        tsteps = res.out_steps
        tvect = [res.steps[kt].time for kt in tsteps]
    except:
        res = zr(pathname,prob)
        res.read_rcf()
        res.read_his()

        tvect = [2,2.99,3,3.99,4,4.99,5,5.99,6,6.99,7,8]
        tsteps = []
        tsteps_plot = []
        for kt,step in enumerate(res.steps):
            if step.time in tvect and step.sf==0 and step.conv_status==-1:
                tsteps.append(kt)
                
        res.out_steps = tsteps
        res.read_dat()
        res.read_s03()
        pickle.dump(res, open(prob+'_etais.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    ylev = set()
    for ke in range(res.nTrusses):
        inel = res.truss.inel[ke]
        y = 0.5*sum([res.coords[1][kn-1] for kn in inel])
        ylev.add(y)
    ylev = sorted(list(ylev),reverse=True)
    
    elists = [[] for k in ylev]
    
    for ke in range(res.nTrusses):
        inel = res.truss.inel[ke]
        y = 0.5*sum([res.coords[1][kn-1] for kn in inel])
        flag = False
        for ky,yy in enumerate(ylev):
            if abs(y-yy)<1e-3:
                elists[ky].append(ke)
                flag = True
        if not flag:
            print('Truss at %1.2f not in ylev-list'%(y))

    N = [[1e10 for ke in range(len(elists[kl]))] for kl in range(len(elists))]
    for kt in tsteps:
        step = res.steps[kt]
        for kl in range(len(elists)):
            for kke,ke in enumerate(elists[kl]):
                N[kl][kke] = min(N[kl][kke],step.truss.force[ke])

    fig = plt.figure(figsize=(20,10))
    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
    fig.text(0.1,0.8,u'$N_{min}$ [kN] über alle Bauphasen',
             size=18,bbox=bbox_props)
    fig.text(0.01,0.01,prob,size=8)

    stit = [u'Spriesslage %i (%1.2f MüM)'%(k+1,ylev[k]) for k in range(len(elists))]

    for kl in range(len(elists)):
        ax = fig.add_subplot(2,3,kl+2)
        ax.axis('off')
        for kke in range(len(elists[kl])):
##            mat = res.truss.mat[elists[kl][kke]]
##            As = 
            inel = res.truss.inel[elists[kl][kke]]
            X = [res.coords[0][kn-1] for kn in inel]
            Z = [-res.coords[2][kn-1] for kn in inel]
            if X[1]-X[0]==0:
                ang = 90
            else:
                ang = math.atan((Z[1]-Z[0])/(X[1]-X[0]))#/math.pi*180
            ax.plot(X,Z,'k')
            ax.annotate('%1.0f'%(N[kl][kke]),xy=(np.mean(X),np.mean(Z)),
                        rotation=0,ha='center',va='center',
                        bbox=dict(boxstyle='round,pad=0.2',fc='w',ec='k'))
####            if elists[kl][kke] in [34,47,8,60,21]:
##            if elists[kl][kke] in [5,18,31,44,57]:
##                ax.annotate('%1.0f kN/m'%(math.cos(ang)*N[kl][kke]/4.39),
##                            xy=(np.mean(X),np.mean(Z)+2),
##                            rotation=0,ha='center',va='center',
##                            bbox=dict(boxstyle='round,pad=0.2',fc='w',ec='k'))
####            elif elists[kl][kke] in [35,48,9,61,22]:
##            elif elists[kl][kke] in [6,19,32,45,58]:
##                ax.annotate('%1.0f kN/m'%(math.cos(ang)*N[kl][kke]/4.02),
##                            xy=(np.mean(X),np.mean(Z)+2),
##                            rotation=0,ha='center',va='center',
##                            bbox=dict(boxstyle='round,pad=0.2',fc='w',ec='k'))
                
        ax.set_title(stit[kl],
                     size=16)
##        ax.set_xlim(-1,91)
##        ax.set_ylim(-70,1)
##        ax.set_xlabel('z [m]',size=16)
##        ax.set_ylabel('x [m]',size=16)
        ax.set_aspect('equal')

    fig.tight_layout()
    fig.savefig(prob+'_Netais')
    plt.close(fig)


        

            




