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


def plot_efforts_etais(prob,res,figsize=(8,10),figpath='.'):
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

    N = [[1e10 for ke in range(len(elists[kl]))] for kl in range(len(elists))]
    P0 = [[-1 for ke in range(len(elists[kl]))] for kl in range(len(elists))]
    for kt in tsteps:
        step = res.steps[kt]
        for kl in range(len(elists)):
            for kke,ke in enumerate(elists[kl]):
                N[kl][kke] = min(N[kl][kke],step.truss.force[ke])
                if res.EF[res.truss.EF[ke]][0]+0.01==step.time:
                    P0[kl][kke] = step.truss.force[ke]

    fig = plt.figure(figsize=figsize)
    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
    fig.text(0.02,0.5,u'$P0$ et $N_{min}$ en kN, $N_{min}$ pour toutes les phases',
             size=18,rotation=90,va='center',bbox=bbox_props)
    fig.text(0.01,0.01,prob,size=8)

    stit = [u'N%i (%1.1f msm)'%(k+1,ylev[k]) for k in range(len(elists))]

    for kl in range(len(elists)):
        ax = fig.add_subplot(3,1,kl+1)
        ax.axis('off')
        for kke in range(len(elists[kl])):
##            mat = res.truss.mat[elists[kl][kke]]
##            As = 
            inel = res.truss.inel[elists[kl][kke]]
            X = [res.coords[0][kn-1] for kn in inel]
            Z = [-res.coords[2][kn-1] for kn in inel]
##            if X[1]-X[0]==0:
##                ang = 90
##            else:
            ang = math.atan((Z[1]-Z[0])/(X[1]-X[0]))/math.pi*180
            ax.plot(X,Z,'k')
            ax.annotate('%s: $P0=$%1.0f, $N_{min}=$%1.0f'%(RORdict[res.truss.mat[elists[kl][kke]]],P0[kl][kke],N[kl][kke]),
                        xy=(np.mean(X),np.mean(Z)),
                        rotation=ang,ha='center',va='top',
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
        ax.set_xlim(104,153)
##        ax.set_ylim(-70,1)
##        ax.set_xlabel('z [m]',size=16)
##        ax.set_ylabel('x [m]',size=16)
        ax.yaxis.set_ticklabels([])
        ax.set_aspect('equal')

    fig.tight_layout()
    fig.savefig(figpath+'/'+prob+'_Netais')
    plt.close(fig)

    return figpath+'/'+prob+'_Netais'


        

            




