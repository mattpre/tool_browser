# -*- coding: cp1252 -*-
import os,math,cmath
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import cPickle as pickle
from matplotlib.patches import Circle
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from zsoil_tools import zsoil_results as zr

pathname = '../..'
pblist = ['M1113_3D_v2','M1113_3D_v2_E60%','M1113_3D_v2_Emin','M1113_3D_v2_EminFlGl50%']


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

        tvect = [0,7,9,20,21]
        tsteps = []
        for kt,step in enumerate(res.steps):
            if step.time in tvect and step.sf==0 and step.conv_status==-1:
                tsteps.append(kt)
##        tsteps = [5,6]
                
        res.out_steps = tsteps
        res.read_dat()
    ##        res.read_s02()
        res.read_s00()
        res.read_s04()
        res.read_s07()
        pickle.dump(res, open(prob+'.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    step = res.steps[tsteps[4]]
    step0 = res.steps[tsteps[3]]
    step = res.steps[tsteps[2]]
    step0 = res.steps[tsteps[1]]
  
    labs = [u'paroi clouée Lac',
            u'paroi clouée Montagne',
            'semelle Lac',
            'semelle Montagne']
    #plane: 3 pts
    planes = [[np.array((-3.03890e+01,-3.48212e+01,0.00000e+00)),
               np.array((-3.01788e+01,-3.37582e+01,0.00000e+00)),
               np.array((-3.03890e+01,-3.48212e+01,4.00000e+00))],
              [np.array((-7.79802e+01,-3.46180e+01,0.00000e+00)),
               np.array((-7.81948e+01,-3.35417e+01,0.00000e+00)),
               np.array((-7.79802e+01,-3.46180e+01,4.00000e+00))],
              [np.array((-3.03890e+01,-3.48212e+01,0.00000e+00)),
               np.array((-3.58754e+01,-3.57870e+01,0.00000e+00)),
               np.array((-3.03890e+01,-3.48212e+01,4.00000e+00))],
              [np.array((-7.79802e+01,-3.46180e+01,0.00000e+00)),
               np.array((-7.25022e+01,-3.55823e+01,0.00000e+00)),
               np.array((-7.79802e+01,-3.46180e+01,4.00000e+00))]]
    cosplane = []
    planes2 = []
    for p in planes:
        v0 = p[1]-p[0]
        v0 /= np.linalg.norm(v0)
        v1 = p[2]-p[0]
        v1 /= np.linalg.norm(v1)
        n = np.cross(v0,v1)
        planes2.append([p[0],v0,v1,n])
    planes2[0][1] *= -1
    planes2[2][1] *= -1
    planes2[0][3] *= -1
    planes2[3][3] *= -1
    elists = [[] for k in labs]
    for ke in range(res.nContacts):
        if res.cnt.type[ke]==1 and res.cnt.EF[ke]==20:
            inel = res.cnt.inel[ke][:4]
            xm = 0.25*sum([res.coords[0][kn-1] for kn in inel])
            ym = 0.25*sum([res.coords[1][kn-1] for kn in inel])
            zm = 0.25*sum([res.coords[2][kn-1] for kn in inel])
            for kp in range(len(planes)):
                p = planes2[kp]
                v = np.array((xm,ym,zm))-p[0]
                if abs(np.dot(v,p[3]))<1e-3:
                    elists[kp].append(ke)
##            if zm<-0.01 and ym<-30:
##                if xm<-77.98:
##                    elists[0].append(ke)
##                elif xm<-72.5 and ym<-34.6:
##                    elists[2].append(ke)
##                elif xm>-35.87 and ym<-34.82:
##                    elists[3].append(ke)
##                elif xm>-30.39:
##                    elists[1].append(ke)
   
    # Sohlpressung:    
    fig = plt.figure(figsize=(16,10))
    bbox_props = dict(boxstyle='round',fc="w", ec="k", lw=2)
    bb2 = dict(boxstyle='square', fc="w",ec='none',alpha=0.5)
    fig.text(0.01,0.01,prob,size=8)
    patches = []
    cvect = []
    minmax = [[[1e10,-1e10],[1e10,-1e10]] for k in range(len(labs))]
    AX = []

    patches = [[] for k in labs]
    cvect = [[] for k in labs]
    unvect = [[] for k in labs]
    subplotorder = [0,3,1,2]
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(3,4,width_ratios=[1,4,4,1])
    for kl in range(len(labs)):
        AX.append([fig.add_subplot(gs[subplotorder[kl]+4*k]) for k in range(3)])
        plane = planes2[kl]
        cosplane.append(plane[3][0])
        for kke,ke in enumerate(elists[kl]):
            xy = np.zeros(4*2)
            xy.shape = (4,2)
            inel = res.cnt.inel[ke][:4]
            un = []
            for kkn,kn in enumerate(inel):
                v = np.array([res.coords[kk][kn-1] for kk in range(3)])
                xy[kkn][0] = np.dot(v-plane[0],plane[1])
                xy[kkn][1] = np.dot(v-plane[0],plane[2])
                minmax[kl][0][0] = min(minmax[kl][0][0],xy[kkn][0])
                minmax[kl][0][1] = max(minmax[kl][0][1],xy[kkn][0])
                minmax[kl][1][0] = min(minmax[kl][1][0],xy[kkn][1])
                minmax[kl][1][1] = max(minmax[kl][1][1],xy[kkn][1])
                uv = np.array([step.nodal.disp[kk][kn-1]
                               -step0.nodal.disp[kk][kn-1] for kk in range(3)])
                un.append(np.dot(uv,plane[3]))
            poly = Polygon(np.array(xy))
            patches[kl].append(poly)
            p_cnt = 0.25*sum([step.cnt.stress[2][ke][kgp]
                              -step0.cnt.stress[2][ke][kgp] for kgp in range(4)])
            cvect[kl].append(p_cnt)
            unvect[kl].append(0.25*sum(un)*1000)
            
##            x = [res.coords[0][kn-1] for kn in inel]
##            z = [res.coords[2][kn-1] for kn in inel]
##            a = (max(x)-min(x))*(max(z)-min(z))
##            A += a
##            for kgp in range(4):
##                sy += (step.cnt.stress[2][ke][kgp]
##                       -step0.cnt.stress[2][ke][kgp])*a*0.25
        
    for kax,axs in enumerate(AX):
        pc = PatchCollection(patches[kax],edgecolors=('none',),
                             antialiased=False)
        pc.set_array(np.array(cvect[kax]))
        axs[0].add_collection(pc)
        axs[0].set_xlim(minmax[kax][0])
        axs[0].set_ylim(minmax[kax][1])
        cb = fig.colorbar(pc,ax=axs[0])
        cb.set_label(u'Pression %s [kPa]'%(labs[kax]))
        cb.ax.tick_params()
        axs[0].set_aspect('equal')

        area = abs((minmax[kax][0][1]-minmax[kax][0][0])
                   *(minmax[kax][1][1]-minmax[kax][1][0]))
        axs[0].annotate('$F_h=$%1.1f kN\n$F_v=$%1.1f kN'%(np.mean(cvect[kax])*area*cosplane[kax],
                                                 np.mean(cvect[kax])*area*(1.-cosplane[kax]**2)),
                        xy=(0.5*sum(minmax[kax][0]),0.5*sum(minmax[kax][1])),
                        ha='center',va='center',rotation=90,size=12,
                        bbox=bb2)
        
        pc = PatchCollection(patches[kax],edgecolors=('none',),
                             antialiased=False)
        pc.set_array(np.array(unvect[kax]))
        axs[1].add_collection(pc)
        axs[1].set_xlim(minmax[kax][0])
        axs[1].set_ylim(minmax[kax][1])
        cb = fig.colorbar(pc,ax=axs[1])
        cb.set_label(u'Dépl. %s [mm]'%(labs[kax]))
        cb.ax.tick_params()
        axs[1].set_aspect('equal')
        
        pc = PatchCollection(patches[kax],edgecolors=('none',),
                             antialiased=False)
        modules = [cvect[kax][kk]/unvect[kax][kk] for kk in range(len(cvect[kax]))]
        pc.set_array(np.array(modules))
        axs[2].add_collection(pc)
        axs[2].set_xlim(minmax[kax][0])
        axs[2].set_ylim(minmax[kax][1])
        cb = fig.colorbar(pc,ax=axs[2])
        cb.set_label(u'$k_B$ %s [MN/m3]'%(labs[kax]))
        cb.ax.tick_params()
        axs[2].set_aspect('equal')

        axs[2].annotate('%1.1f MN/m3'%(sum(modules)/len(modules)),
                        xy=(0.5*sum(minmax[kax][0]),0.5*sum(minmax[kax][1])),
                        ha='center',va='center',rotation=90,size=12,
                        bbox=bb2)
        
    fig.tight_layout()
    tstr = 'T%1.1f-T%1.1f'%(step.time,step0.time)
    fig.text(0.2,0.98,tstr,size=14,
             ha='left',va='top',bbox=bbox_props)
    import re
    a = re.sub('\.','_',tstr)
    fig.savefig(prob+'_ressorts_%s'%(a))
    plt.close()
   








