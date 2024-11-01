# -*- coding: cp1252 -*-
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pickle
from matplotlib.patches import Circle
from matplotlib.patches import Polygon
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from scipy import interpolate
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from zsoil_tools import zsoil_results as zr
from zsoil_tools import postpro_lib as pl

def get_quad_area_3d(xx,yy,zz):
    p = [np.array([xx[k],yy[k],zz[k]]) for k in range(4)]
    # tri 0:
    v0 = p[1]-p[0]
    v1 = p[3]-p[0]
    a0 = 0.5*np.linalg.norm(np.cross(v0,v1))
    # tri 1:
    v0 = p[2]-p[1]
    v1 = p[3]-p[1]
    a1 = 0.5*np.linalg.norm(np.cross(v0,v1))
    return a0+a1

pathname = '..'
pblist = ['M1172_vol_v3_2(noPrec)',
          'M1172_vol_v3_2(noPrec)_pess',
          'M1172_vol_v3_2(noPrec)_piles16m']
pblist = ['M1172_vol_v3_1(chariot)',
          'M1172_vol_v3_1(chariot)_eau2']
plabs = ['Base',
         'eau2']

# axe neutre:
y0 = 0.22288
coupes = [(u'mi-travée',19.55),
          (u'appui',-1.45)]
disp_nodes = [(u'mi-travée',np.array([19.55,8.91,0])),
              (u'appui',np.array([-1.45,8.91,0]))]

TX = [[]for _ in pblist]
MM = [[[] for _ in coupes] for _ in pblist]
TT = [[[] for _ in coupes] for _ in pblist]
DY = [[[] for _ in coupes] for _ in pblist]
M0 = [[[[] for _ in range(5)] for _ in coupes] for _ in pblist]
for kf,prob in enumerate(pblist):
    figname = prob
    try:
        res = pickle.load(open(prob+'_tablier.p', "rb" ))
        print(prob+' loaded')
        tsteps = res.out_steps
        tvect = [res.steps[kt].time for kt in tsteps]
    except:
        res = zr(pathname,prob)
        res.read_rcf()
        res.read_his()
        res.ref_vect2D = [1,0,0]

        tvect = [5,6,6.5,7,8,9,9.25,9.5,9.75,10]
        tsteps = []
        for kt,step in enumerate(res.steps):
            if step.time>3.9 and step.sf==0 and step.conv_status==-1:
##            if step.time in tvect and step.sf==0 and step.conv_status==-1:
                tsteps.append(kt)
##        tsteps = [5,6]
                
        res.out_steps = tsteps
        res.read_dat()
        res.read_s02()
        res.read_s00()
        res.read_s04()
        pickle.dump(res, open(prob+'_tablier.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)


    cul_nodes = set()
    quads = []
    for ke in range(res.nVolumics):
        if res.vol.mat[ke] == 14:
            inel = res.vol.inel[ke]
            xx = [res.coords[0][kn-1] for kn in inel]
            yy = [res.coords[1][kn-1] for kn in inel]
            zz = [res.coords[2][kn-1] for kn in inel]
            if max(zz)>-0.01:
                quad = []
                for kk in range(8):
                    if zz[kk]>-0.01:
                        quad.append([xx[kk],yy[kk]])
                center = np.mean(quad,0)
                quad = [(q,np.angle(q[0]-center[0]+1j*(q[1]-center[1]))) for q in quad]
                temp = sorted(quad,key=lambda v:v[1])
                quads.append([t[0] for t in temp])
            for kn in inel:
                crd = [res.coords[kk][kn-1] for kk in range(3)]
                if not (abs(crd[0]+1.45)<1e-3 and abs(crd[1]-7.75)<1e-3):
##                if not (abs(res.coords[0][kn-1]+1.45)<1e-3 and res.coords[1][kn-1]>7.74):
##                    if not (res.coords[0][kn-1]<-5.8 and res.coords[1][kn-1]>8.0):
                    cul_nodes.add(kn-1)
    cul_nodes = list(cul_nodes)
    

    tpairs = [(4,6),(6,7),(7,9),(9,10)]
    tlabs = [u'Poids de la culée supérieure + pont (avant clavage)',
             u'Précontrainte',
             u'Remblayage + dalle de transition',
             u'Charges de trafic']
    maxval = [[[1e10,-1e10],[1e10,-1e10],0] for _ in tpairs]
    tvect = [res.steps[kt].time for kt in tsteps]

    of = open(prob+'_disp.csv','w')
    of.write(';')
    for kk,tpair in enumerate(tpairs):
        step0 = res.steps[tsteps[tvect.index(tpair[0])]]
        step1 = res.steps[tsteps[tvect.index(tpair[1])]]
        of.write(';T%1.0f-T%1.0f'%(step1.time,step0.time))

        X = []
        Y = []
        U = []
        V = []
        for kn in cul_nodes:
            maxval[kk][0][0] = min(maxval[kk][0][0],step1.nodal.disp[0][kn]-step0.nodal.disp[0][kn])
            maxval[kk][0][1] = max(maxval[kk][0][1],step1.nodal.disp[0][kn]-step0.nodal.disp[0][kn])
            maxval[kk][1][0] = min(maxval[kk][1][0],step1.nodal.disp[1][kn]-step0.nodal.disp[1][kn])
            maxval[kk][1][1] = max(maxval[kk][1][1],step1.nodal.disp[1][kn]-step0.nodal.disp[1][kn])
            maxval[kk][2] = max(maxval[kk][2],np.linalg.norm(np.array([step1.nodal.disp[0][kn]-step0.nodal.disp[0][kn],
                                                                       step1.nodal.disp[1][kn]-step0.nodal.disp[1][kn]])))
            crd = [res.coords[kk][kn] for kk in range(3)]
            if crd[2]>-0.01:
                X.append(crd[0])
                Y.append(crd[1])
                U.append(step1.nodal.disp[0][kn]-step0.nodal.disp[0][kn])
                V.append(step1.nodal.disp[1][kn]-step0.nodal.disp[1][kn])


        fig = plt.figure(figsize=(6,8))
        bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
        fig.text(0.01,0.01,prob,size=8)

        ax = fig.add_subplot(111)

        for quad in quads:
            ax.add_patch(Polygon(np.array(quad),alpha=1,ec='none',fc='#ffff80',aa='none',zorder=0))

        cell_text = [['MIN-UX','[mm]','%1.1f'%(maxval[kk][0][0]*1e3)],
                     ['MAX-UX','[mm]','%1.1f'%(maxval[kk][0][1]*1e3)],
                     ['MIN-UY','[mm]','%1.1f'%(maxval[kk][1][0]*1e3)],
                     ['MAX-UY','[mm]','%1.1f'%(maxval[kk][1][1]*1e3)],
                     ['MAX |U|','[mm]','%1.1f'%(maxval[kk][2]*1e3)]]
        table = plt.table(cellText=cell_text,
                          bbox=[0.6,0.4,0.3,0.2])
        table.auto_set_column_width(col=[0,1,2])
        table.set_fontsize(10)
        
        fig.text(0.5,0.04,'T=%1.0f-T=%1.0f: %s'%(step1.time,step0.time,tlabs[kk]),size=12,ha='center')

        ax.quiver(X,Y,U,V,color='#0080c0',scale=0.05)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_ylim([-0.5,9])

        fig.tight_layout()
        fig.savefig(prob+'_arrows_T%1.0f-T%1.0f'%(step1.time,step0.time))
        
    of.write('\nMIN-UX;[mm]')
    for kk,tpair in enumerate(tpairs):
        of.write(';%1.6e'%(maxval[kk][0][0]))
    of.write('\nMAX-UX;[mm]')
    for kk,tpair in enumerate(tpairs):
        of.write(';%1.6e'%(maxval[kk][0][1]))
    of.write('\nMIN-UY;[mm]')
    for kk,tpair in enumerate(tpairs):
        of.write(';%1.6e'%(maxval[kk][1][0]))
    of.write('\nMAX-UY;[mm]')
    for kk,tpair in enumerate(tpairs):
        of.write(';%1.6e'%(maxval[kk][1][1]))
    of.write('\nMAX |U|;[mm]')
    for kk,tpair in enumerate(tpairs):
        of.write(';%1.6e'%(maxval[kk][2]))
    of.close()

   
