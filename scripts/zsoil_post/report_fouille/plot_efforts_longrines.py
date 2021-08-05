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
import pickle
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator

HEBdict = {10:'2xHEB400',
           11:'2xHEB600'}


def plot_efforts_longrines(prob,res,figsize=(8,8),figpath='.'):
    tsteps = res.out_steps
    
    ylev = [703.5,699]#,694.5]
    
    elists = [[] for k in ylev]
    
    for ke in range(res.nBeams):
        inel = res.beam.inel[ke]
        y = 0.5*sum([res.coords[1][kn-1] for kn in inel])
        flag = False
        for ky,yy in enumerate(ylev):
            if abs(y-yy)<1e-3:
                elists[ky].append(ke)
                flag = True

    N = [[[1e10,-1e10] for ke in range(len(elists[kl]))] for kl in range(len(elists))]
    M = [[[[1e10,1e10],[-1e10,-1e10]] for ke in range(len(elists[kl]))] for kl in range(len(elists))]
    Q = [[[1e10,-1e10] for ke in range(len(elists[kl]))] for kl in range(len(elists))]
    scn = 0.002
    scm = -0.008
    scq = 0.01
    for kt in tsteps:
        step = res.steps[kt]
        for kl in range(len(elists)):
            for kke,ke in enumerate(elists[kl]):
                inel = res.beam.inel[ke]
                X = [res.coords[0][kn-1] for kn in inel]
                Z = [-res.coords[2][kn-1] for kn in inel]
                dl = 0.5*((X[1]-X[0])**2+(Z[1]-Z[0])**2)**0.5
                N[kl][kke][0] = min(N[kl][kke][0],step.beam.force[0][ke])
                N[kl][kke][1] = max(N[kl][kke][1],step.beam.force[0][ke])
                M[kl][kke][0][0] = min(M[kl][kke][0][0],step.beam.moment[1][ke]-dl*step.beam.force[2][ke])
                M[kl][kke][0][1] = min(M[kl][kke][0][1],step.beam.moment[1][ke]+dl*step.beam.force[2][ke])
                M[kl][kke][1][0] = max(M[kl][kke][1][0],step.beam.moment[1][ke]-dl*step.beam.force[2][ke])
                M[kl][kke][1][1] = max(M[kl][kke][1][1],step.beam.moment[1][ke]+dl*step.beam.force[2][ke])
                Q[kl][kke][0] = min(Q[kl][kke][0],step.beam.force[2][ke])
                Q[kl][kke][1] = max(Q[kl][kke][1],step.beam.force[2][ke])
    minmax = [[min([min([min(N[kl][kke]) for kke in range(len(elists[kl]))]) for kl in range(len(elists))]),
               max([max([max(N[kl][kke]) for kke in range(len(elists[kl]))]) for kl in range(len(elists))])],
              [min([min([min(min(M[kl][kke][0]),min(M[kl][kke][1])) for kke in range(len(elists[kl]))]) for kl in range(len(elists))]),
               max([max([max(max(M[kl][kke][0]),max(M[kl][kke][1])) for kke in range(len(elists[kl]))]) for kl in range(len(elists))])],
              [min([min([min(Q[kl][kke]) for kke in range(len(elists[kl]))]) for kl in range(len(elists))]),
               max([max([max(Q[kl][kke]) for kke in range(len(elists[kl]))]) for kl in range(len(elists))])]]

    FIG = [plt.figure(figsize=figsize) for kk in range(3)]
    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
    AX = []
    norms = []
    ticks = []
    dround = [500,100,100]
    for kf,fig in enumerate(FIG):
        fig.text(0.01,0.01,prob,size=8)
        ax = fig.add_subplot(111)
        ax.axis('off')
        AX.append(ax)
        ticks.append(MultipleLocator(dround[kf]).tick_values(minmax[kf][0],minmax[kf][1]))
        norms.append(colors.BoundaryNorm(ticks[kf],len(ticks[kf])))

    labs = ['Effort normal','Moment','Effort tranchant']
    labs1 = ['N','M','Q']
    units = ['kN','kNm','kN']
    stit = [u'%s, N%i (%1.1f msm)'%(labs[kl],k+1,ylev[k]) for k in range(len(elists))]

    pcs = [0,0,0]
    for kl in range(len(elists)):
        patches = [[],[],[]]
        cvect = [[],[],[]]
        mvect = []
        for kke in range(len(elists[kl])):
##            mat = res.truss.mat[elists[kl][kke]]
##            As = 
            inel = res.beam.inel[elists[kl][kke]]
            X = [np.array([res.coords[0][kn-1],-res.coords[2][kn-1]]) for kn in inel]
            n = X[1]-X[0]
            n = np.array([-n[1],n[0]])
            n /= np.linalg.norm(n)
            X += np.array([0,-20*kl])
##            if X[1]-X[0]==0:
##                ang = 90
##            else:
##            ang = math.atan((Z[1]-Z[0])/(X[1]-X[0]))/math.pi*180
            # effort normal:
            AX[0].plot([X[0][0],X[1][0]],[X[0][1],X[1][1]],'k')
            XX = [X[0],X[0]+n*scn*N[kl][kke][1],
                  X[1]+n*scn*N[kl][kke][1],X[1]]
            poly = Polygon(XX)
            patches[0].append(poly)
            cvect[0].append(N[kl][kke][1])
            XX = [X[0],X[0]+n*scn*N[kl][kke][0],
                  X[1]+n*scn*N[kl][kke][0],X[1]]
            poly = Polygon(XX)
            patches[0].append(poly)
            cvect[0].append(N[kl][kke][0])
            # moment:
            AX[1].plot([X[0][0],X[1][0]],[X[0][1],X[1][1]],'k')
            XX = [X[0],X[0]+n*scm*M[kl][kke][1][0],
                  X[1]+n*scm*M[kl][kke][1][1],X[1]]
            poly = Polygon(XX)
            patches[1].append(poly)
            cvect[1].append(np.mean(M[kl][kke][1]))
            mvect.append(M[kl][kke][1][0])
            mvect.append(M[kl][kke][1][1])
            XX = [X[0],X[0]+n*scm*M[kl][kke][0][0],
                  X[1]+n*scm*M[kl][kke][0][1],X[1]]
            poly = Polygon(XX)
            patches[1].append(poly)
            cvect[1].append(np.mean(M[kl][kke][0]))
            mvect.append(M[kl][kke][0][0])
            mvect.append(M[kl][kke][0][1])
            # effort tranchant:
            AX[2].plot([X[0][0],X[1][0]],[X[0][1],X[1][1]],'k')
            XX = [X[0],X[0]+n*scq*Q[kl][kke][1],
                  X[1]+n*scq*Q[kl][kke][1],X[1]]
            poly = Polygon(XX)
            patches[2].append(poly)
            cvect[2].append(Q[kl][kke][1])
            XX = [X[0],X[0]+n*scq*Q[kl][kke][0],
                  X[1]+n*scq*Q[kl][kke][0],X[1]]
            poly = Polygon(XX)
            patches[2].append(poly)
            cvect[2].append(Q[kl][kke][0])
        for kk in range(3):
            vmin = min(cvect[kk])
            ind = cvect[kk].index(vmin)
            if kk==1:
                AX[kk].annotate('%1.0f %s'%(min(mvect),units[kk]),
                                xy=(np.mean(patches[kk][ind].xy[1:3],0)),
                                ha='center',va='center')
            else:
                AX[kk].annotate('%1.0f %s'%(vmin,units[kk]),
                                xy=(np.mean(patches[kk][ind].xy[1:3],0)),
                                ha='center',va='center')
            vmax = max(cvect[kk])
            ind = cvect[kk].index(vmax)
            if kk==1:
                AX[kk].annotate('%1.0f %s'%(max(mvect),units[kk]),
                                xy=(np.mean(patches[kk][ind].xy[1:3],0)),
                                ha='center',va='center')
            else:
                AX[kk].annotate('%1.0f %s'%(vmax,units[kk]),
                                xy=(np.mean(patches[kk][ind].xy[1:3],0)),
                                ha='center',va='center')
            pc = PatchCollection(patches[kk],edgecolors='none',cmap=plt.cm.get_cmap('viridis'))
            pc.set_array(np.array(cvect[kk]))
            pc.set_clim(vmin=minmax[kk][0],vmax=minmax[kk][1])
            pcs[kk] = pc
            AX[kk].add_collection(pc)
            AX[kk].annotate(u'N%i (%1.1f msm)'%(kl+1,ylev[kl]),xy=(128,-111-20*kl),ha='center',
                            size=14)
    for kk in range(3):
        cb = FIG[kk].colorbar(pcs[kk],norm=norms[kk],ticks=ticks[kk])
        cb.set_label('%s [%s]'%(labs[kk],units[kk]),size=16)
        cb.ax.tick_params(labelsize=14)
    for ax in AX:
        ax.set_xlim(100,157)
        ax.set_ylim(-175,-105)
        ax.yaxis.set_ticklabels([])
        ax.set_aspect('equal')

    figpaths = []
    for kk,fig in enumerate(FIG):
        fig.tight_layout()
        fig.savefig(figpath+'/'+prob+'_Longrines_%s'%(labs1[kk]))
        plt.close(fig)
        figpaths.append(figpath+'/'+prob+'_Longrines_%s'%(labs1[kk]))

    return figpaths
        

            




