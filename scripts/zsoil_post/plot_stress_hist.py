# -*- coding: cp1252 -*-
# @description Extracting time series from zsoil results
# @input zsoil results
# @output time history in txt-file
# @author Matthias Preisig
# @date 2019/01/11
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

from zsoil_tools import zsoil_results as zr

pathname = '..'

pblist = ['Oedo_Methode_G2_16.02_16.02_OCR_G2']

for kf,prob in enumerate(pblist):
    res = zr(pathname,prob)
    res.read_rcf()
    res.read_his()
    res.read_eda()
    tsteps = range(0,len(res.steps))
    tsteps = []
    for kt,step in enumerate(res.steps):
        if step.conv_status==-1:
            tsteps.append(kt)
    res.out_steps = tsteps
    res.read_dat()
    res.read_LTF()
    res.read_s01('/v500')

    res.compute_invariants()
    
    fig = plt.figure(figsize=(14,10))
    ax = fig.add_subplot(111)

    for ke in range(0,res.nElements,1):
        matnum = res.vol.mat[ke]
        if matnum==2:
            for km,mat in enumerate(res.materials):
                if mat.number==matnum:
                    break
            phi = res.materials[km].nonlinear['Friction angle']
            c = res.materials[km].nonlinear["Effective cohesion c'"]
            Mp = 6*math.sin(phi/180*math.pi)/(3-math.sin(phi/180*math.pi))

            p = []
            q = []
            for kt in tsteps:
                step = res.steps[kt]
                p.append(-step.vol.invar[2][ke])
                q.append(step.vol.invar[3][ke])

            ax.plot([0,1.5*p[-1]],[0,1.5*p[-1]*Mp])

            points = np.array([p,q]).T.reshape(-1,1,2)
            segments = np.concatenate([points[:-1],points[1:]],axis=1)
            norm = plt.Normalize(0,len(p)-1)
            lc = LineCollection(segments,norm=norm)
            lc.set_array(np.array(range(0,len(p))))
##            lc.set_linewidth(2)
            line = ax.add_collection(lc)
##        ax.plot(p,q,'.')

    ax.set_xlabel(u'p [kPa]')
    ax.set_ylabel(u'q [kPa]')
##    ax.legend(prop={'size':'small'},loc=2)
    ax.grid(b=True,which='both',axis='x')
    ax.grid(b=True,which='major',axis='y')

    fig.tight_layout()
    fig.savefig(prob+'_p-q.png')



