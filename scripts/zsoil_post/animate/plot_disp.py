# -*- coding: cp1252 -*-
# @description Plotting deformed shape
# @input zsoil results
# @output png figure
# @author Matthias Preisig
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import re


from zsoil_tools import zsoil_results as zr
import time




pathname = '..'
pblist = ['bloc_v2']

t0 = 2
for kf,prob in enumerate(pblist):
    res = zr(pathname,prob)
    res.read_rcf()
    res.read_his()
    tsteps = []
    for kt,step in enumerate(res.steps):
        if step.conv_status in [-1]:
            if abs(float('%1.0f'%(step.time*20))-step.time*20)<1e-3:
                tsteps.append(kt)
            
    res.out_steps = tsteps
    res.read_dat()
    res.read_s00('/v500')
##    res.read_s00(res_type='reactions')
##    res.read_s01('/v500')
    

    inel = np.zeros(res.nVolumics*2*3)
    inel.shape = (res.nVolumics*2,3)
    kv0 = [0,1,2]
    kv1 = [0,2,3]
    for ke in range(res.nVolumics):
        if res.vol.EF[ke] in [0,1]:
            for kk in range(3):
                inel[ke*2][kk] = res.vol.inel[ke][kv0[kk]]-1
                inel[ke*2+1][kk] = res.vol.inel[ke][kv1[kk]]-1

    for kn in range(res.nNodes):
        if abs(res.coords[0][kn]-1.24317e+03)<0.1 and res.coords[1][kn]<965.1:
            kn0 = kn
            break

    cmap = plt.cm.jet
    maxval = 0.01#max(val)
    norm = mpl.colors.Normalize(vmin=0,vmax=maxval)

    for kkt in range(len(tsteps)):

        fig = plt.figure(figsize=(6,6))
        bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
        fig.text(0.01,0.01,prob,size=8)
        ax = fig.add_subplot(111)
        
        step = res.steps[tsteps[kkt]]

        bounds = [[1e10,-1e10],
                  [1e10,-1e10]]
    
        scale = 1
        dx = [step.nodal.disp[0][kn] for kn in range(res.nNodes)]
        dy = [step.nodal.disp[1][kn] for kn in range(res.nNodes)]
        x = np.array([res.coords[0][kn]+dx[kn]*scale for kn in range(res.nNodes)])
        y = np.array([res.coords[1][kn]+dy[kn]*scale for kn in range(res.nNodes)])
        triang = tri.Triangulation(x,y,inel)

        val = []
        for kn in range(res.nNodes):
            val.append(math.sqrt(dx[kn]**2+
                                 dy[kn]**2))
##        norm = mpl.colors.Normalize(vmin=0,vmax=max(val))
        a = ax.tripcolor(triang,val,cmap=cmap,norm=norm,shading='gouraud')
        ax.annotate('%1.3f m'%(max(val)),xy=(-190,1210))
        
        ax.set_xlim(-2,10)
        ax.set_ylim(0,12)
##        ax.set_ylabel(u'[m ü. M.]')
##        ax.set_ylim([410,510])
        ax.grid('on')
        ax.set_aspect('equal')

        fig.tight_layout()
        fig.savefig(prob+re.sub('\.','_','_%1.2f'%(step.time))+'.jpeg')
        plt.close(fig)



