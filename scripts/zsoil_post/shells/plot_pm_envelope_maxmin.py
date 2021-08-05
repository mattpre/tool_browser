# -*- coding: cp1252 -*-
# @description Plotting envelopes of shell forces for vertical walls
# @input zsoil results
# @output png
# @author Matthias Preisig
# @date 2017/11/08
import os,math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cPickle as pickle
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from zsoil_tools import zsoil_results as zr
import time

#############################################################
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator
def plot_cmap(ax,cvect,patches,ncol,minmax=0,labstr=0):
    if minmax==0:
        minima = min(cvect)
        maxima = max(cvect)
    else:
        minima = minmax[0]
        maxima = minmax[1]
    dval = (maxima-minima)/(ncol-2)
    nice = math.pow(10,math.floor(math.log10(dval)))
    if dval/nice<1.4:
        dround = nice
    elif dval/nice<3:
        dround = nice*2
    elif dval/nice<7:
        dround = nice*5
    else:
        dround = nice*10
        
    ticks = MultipleLocator(dround).tick_values(minima,maxima)

    cmap = plt.cm.get_cmap('jet')
    norm = colors.BoundaryNorm(ticks,len(ticks))
    pc = PatchCollection(patches,edgecolors=('none',),cmap=cmap)
    pc.set_clim(vmin=minima,vmax=maxima)
    pc.set_array(np.array(cvect))
    ax.add_collection(pc)
    cb = fig.colorbar(pc,cmap=cmap,norm=norm,boundaries=ticks)
    cb.set_ticks(ticks)
    cb.ax.tick_params(labelsize=16)
    if not labstr==0:
        cb.set_label(labstr,size=16)
#############################################################
        

pathname = '//192.168.1.51/Mandats sur H RAID0/M1071_HallesMorges'
pblist = ['M1071_MorgesHalles_3D_TN378_shellHinge_fiche_phase2_DoF_+E_p0_etaiCoin_radier_beamLoads_cages_colonnes_geolOct17_etayv3Forasol']
pblist = ['M1071_MorgesHalles_3D_TN378_shellHinge_fiche_phase2_DoF_+E_p0_etaiCoin_radier_beamLoads_cages_colonnes_geolOct17_NP_etayv3Forasol']

for kf,prob in enumerate(pblist):
    figname = prob
    try:
        res = pickle.load(open(prob+'.p', "rb" ))
        print prob+' loaded'
        tsteps = res.out_steps
        tvect = [res.steps[kt].time for kt in tsteps]
    except:
        res = zr(pathname,prob)
        res.read_rcf()
        res.read_his()

        tvect = [4,5,6,7,8,12]
        tsteps = []
        for kt,step in enumerate(res.steps):
            if step.time in tvect and step.sf==0 and step.conv_status==-1:
                tsteps.append(kt)
##        tsteps = [5,6]
                
        res.out_steps = tsteps
        res.read_dat()
        res.read_s02()
        res.read_s00()
        pickle.dump(res, open(prob+'.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    step0 = res.steps[tsteps[0]]
    step = res.steps[tsteps[1]]

 

    xzpos = []
    xz = []
    elists = []
    M = [[],[]]
    N = [[],[]]
    T = [[],[]]
    for ke in range(res.nShells):
        if res.shell.mat[ke]==21:
            inel = res.shell.inel[ke]
            x = [res.coords[0][kn-1] for kn in inel]
            y = [res.coords[1][kn-1] for kn in inel]
            z = [-res.coords[2][kn-1] for kn in inel]
            xm = 0.25*sum(x)
            ym = 0.25*sum(y)
            zm = 0.25*sum(z)
            posint = int(xm*1000)+int(zm*10)
            if not posint in xzpos:
                xzpos.append(posint)
                if abs(x[0]-x[1])<1e-3 and abs(z[0]-z[1])<1e-3:
                    xz.append([[x[3],x[0]],[z[3],z[0]]])
                else:
                    xz.append([[x[0],x[1]],[z[0],z[1]]])
                elists.append([ke])
                M[0].append(1e10)
                M[1].append(-1e10)
                N[0].append(1e10)
                N[1].append(-1e10)
                T[0].append(1e10)
                T[1].append(-1e10)
            else:
                elists[xzpos.index(posint)].append(ke)

    for kt in tsteps:
        step = res.steps[kt]
        for ke in range(res.nShells):
            for kl,elist in enumerate(elists):
                if ke in elist:
                    break
            if res.shell.mat[ke]==21 :
                M[0][kl] = min(M[0][kl],step.shell.smoment[0][ke])
                M[1][kl] = max(M[1][kl],step.shell.smoment[0][ke])
                N[0][kl] = min(N[0][kl],step.shell.smforce[0][ke])
                N[1][kl] = max(N[1][kl],step.shell.smforce[0][ke])
                T[0][kl] = min(T[0][kl],step.shell.sqforce[0][ke])
                T[1][kl] = max(T[1][kl],step.shell.sqforce[0][ke])

    for kplot in range(3):
        fig = plt.figure(figsize=(15,8))
        bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
        fig.text(0.01,0.01,prob,size=8)

        ax = fig.add_subplot(111)
        if kplot==0:
            sc = 0.01
            val = M
            lab = 'M'
            label = u'Moments usuels (enveloppe sur durée des \ntravaux et hauteur de la paroi) [kNm/m]'
            units = 'kNm/m'
        elif kplot==1:
            sc = 0.01
            val = N
            lab = 'N'
            label = u'Efforts normaux (enveloppe sur durée des \ntravaux et hauteur de la paroi) [kN/m]'
            units = 'kN/m'
        else:
            sc = 0.01
            val = T
            lab = 'T'
            label = u'Efforts tranchants usuels (enveloppe sur durée des \ntravaux et hauteur de la paroi) [kN/m]'
            units = 'kN/m'
        patches = []
        cvect = []
        for kl,elist in enumerate(elists):
            vec = np.array([xz[kl][0][1]-xz[kl][0][0],xz[kl][1][1]-xz[kl][1][0]])
            vec /= np.linalg.norm(vec)
            ax.plot(xz[kl][0],xz[kl][1],'k',linewidth=1)

            poly = Polygon(np.array([[xz[kl][0][0],
                                      xz[kl][0][0]-vec[1]*sc*val[0][kl],
                                      xz[kl][0][1]-vec[1]*sc*val[0][kl],
                                      xz[kl][0][1]],
                                     [xz[kl][1][0],
                                      xz[kl][1][0]+vec[0]*sc*val[0][kl],
                                      xz[kl][1][1]+vec[0]*sc*val[0][kl],
                                      xz[kl][1][1]]]).T)
            patches.append(poly)
            cvect.append(val[0][kl])
            poly = Polygon(np.array([[xz[kl][0][0],
                                      xz[kl][0][0]-vec[1]*sc*val[1][kl],
                                      xz[kl][0][1]-vec[1]*sc*val[1][kl],
                                      xz[kl][0][1]],
                                     [xz[kl][1][0],
                                      xz[kl][1][0]+vec[0]*sc*val[1][kl],
                                      xz[kl][1][1]+vec[0]*sc*val[1][kl],
                                      xz[kl][1][1]]]).T)
            patches.append(poly)
            cvect.append(val[1][kl])
    ##        ax.plot([xz[kl][0][0]-vec[1]*scM*M[0][kl],
    ##                 xz[kl][0][0]-vec[1]*scM*M[1][kl],
    ##                 xz[kl][0][1]-vec[1]*scM*M[1][kl],
    ##                 xz[kl][0][1]-vec[1]*scM*M[0][kl],
    ##                 xz[kl][0][0]-vec[1]*scM*M[0][kl]],
    ##                [xz[kl][1][0]+vec[0]*scM*M[0][kl],
    ##                 xz[kl][1][0]+vec[0]*scM*M[1][kl],
    ##                 xz[kl][1][1]+vec[0]*scM*M[1][kl],
    ##                 xz[kl][1][1]+vec[0]*scM*M[0][kl],
    ##                 xz[kl][1][0]+vec[0]*scM*M[0][kl]],'k',linewidth=0.5)
        plot_cmap(ax,cvect,patches,22,minmax=(min(val[0]),max(val[1])),
                  labstr=label)

        ax.annotate('%smin = %1.0f %s\n%smax = %1.0f %s'%(lab,min(val[0]),units,
                                                          lab,max(val[1]),units),
                    xy=(0.5,0.5),xycoords='axes fraction',ha='center',va='center',size=14)

        ax.set_xlim([-10,140])
        ax.set_ylim([-10,70])
        ax.set_aspect('equal')

        fig.tight_layout()
        fig.savefig(prob+'_pm_'+lab)
   








