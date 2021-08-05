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


RORdict = {15:'609x12',
           16:'711x12.5',
           17:'813x14',
           18:'914x16'}

def plot_efforts_evol_etais(prob,res,figsize=(8,10),figpath='.'):
    tsteps = res.out_steps

    ylev = [703.5,699,694.5]
    
    elists = [[] for k in ylev]
    
    for ke in range(res.nTrusses):
        inel = res.truss.inel[ke]
        y = 0.5*sum([res.coords[1][kn-1] for kn in inel])
        flag = False
        for ky,yy in enumerate(ylev):
            if abs(y-yy)<1e-3:
                elists[ky].append(ke)
                flag = True
##        if not flag:
##            print('Truss at %1.2f not in ylev-list'%(y))

    N = [[[] for ke in range(len(elists[kl]))] for kl in range(len(elists))]
    P0 = [[-1 for ke in range(len(elists[kl]))] for kl in range(len(elists))]
    tx = []
    for kt in tsteps:
        step = res.steps[kt]
        tx.append(step.time)
        for kl in range(len(elists)):
            for kke,ke in enumerate(elists[kl]):
                N[kl][kke].append(step.truss.force[ke])
                if res.EF[res.truss.EF[ke]][0]+0.01==step.time:
                    P0[kl][kke] = step.truss.force[ke]

    fig = plt.figure(figsize=figsize)
    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
##    fig.text(0.02,0.5,u'$P0$ et $N_{min}$ en kN, $N_{min}$ pour toutes les phases',
##             size=18,rotation=90,va='center',bbox=bbox_props)
    fig.text(0.01,0.01,prob,size=8)

    stit = [u'N%i (%1.1f msm)'%(k+1,ylev[k]) for k in range(len(elists))]

    lab = [u'1er (plus court)',u'2ème',u'3ème',u'4ème (plus long)']
    for kl in range(len(elists)):
        ax = fig.add_subplot(3,1,kl+1)
##        ax.axis('off')
        for kke in range(len(elists[kl])):
            kt0 = tx.index(res.EF[res.truss.EF[elists[kl][kke]]][0]+0.01)
            ax.plot(tx[kt0:],N[kl][kke][kt0:],label=lab[kke])
                
        ax.set_title(stit[kl],
                     size=16)
        ax.set_xlim(3,res.steps[tsteps[-1]].time)
##        ax.set_ylim(-70,1)
        ax.set_xlabel('Etapes T',size=16)
        ax.set_ylabel('Effort normal [kN]',size=16)
        ax.grid('on')
        ax.legend()
        ##        ax.set_aspect('equal')

    fig.tight_layout()
    fig.savefig(figpath+'/'+prob+'_Netais_evol')
    plt.close(fig)

    return figpath+'/'+prob+'_Netais_evol'
        

            




