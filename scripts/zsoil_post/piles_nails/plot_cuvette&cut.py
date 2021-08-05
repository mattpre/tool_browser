# -*- coding: cp1252 -*-
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pickle
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.gridspec import GridSpec
colvect = ['#1f77b4','#ff7f0e','#2ca02c']

from zsoil_tools import zsoil_results as zr
from zsoil_tools import postpro_lib as pl

pathname = '..'
pblist = ['M1279_Scheibe_DN600_v2_2',
          # 'M1279_Scheibe_DN600_sigRef_v2_1',
           'M1279_Scheibe_DN600_Emin_v2_2',
           'M1279_Scheibe_DN600_Emax_v2_2',
           'M1279_Scheibe_DN600_fins30m_v2_2',
           'M1279_Scheibe_DN600_DN450axe2_v2_2',
           'M1279_Scheibe_DN600_L-5_v2_2',
          ]
pblist = ['M1279_Scheibe_DN600_fins30m_pieux40m_v2_2']


for kf,prob in enumerate(pblist):
    figname = prob
    try:
        res = pickle.load(open(prob+'.p', "rb" ))
        print(prob+' loaded')
        tsteps = res.out_steps
        tvect = [res.steps[kt].time for kt in tsteps]
    except:
        res = zr(pathname,prob)
        res.read_rcf()
        res.read_his()
    
        tvect = [2,3]
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
        # pickle.dump(res, open(prob+'.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    volnds = set()
    for ke in range(res.nVolumics):
        inel = res.vol.inel[ke]
        for kn in inel:
            volnds.add(kn-1)
    volnds = list(volnds)
    zyndlist = []
    zyposlist = []
    zyposlist1 = []
    for kn in volnds:
        pos = (-res.coords[2][kn],res.coords[1][kn])
        pos1 = res.coords[2][kn]*1e6+res.coords[1][kn]
        if pos1 in zyposlist1:
            ind = zyposlist1.index(pos1)
            zyndlist[ind].append(kn)
        else:
            ind = len(zyposlist)
            zyposlist.append(pos)
            zyposlist1.append(pos1)
            zyndlist.append([kn])

    triang = mtri.Triangulation([v[0] for v in zyposlist],
                                [v[1] for v in zyposlist])

    step = res.steps[tsteps[1]]
    step0 = res.steps[tsteps[0]]

    vals = []
    for nlist in zyndlist:
        vals.append(np.mean([(step.nodal.disp[1][kn]-step0.nodal.disp[1][kn])*1000 for kn in nlist]))

    fig = plt.figure(figsize=(20,6))
    fig.text(0.01,0.01,prob,size=8)
    gs = GridSpec(2,2,width_ratios=[0.96,0.04],height_ratios=[1,3])
    ax0 = fig.add_subplot(gs[0,0])
    ax1 = fig.add_subplot(gs[1,0])
    cax = fig.add_subplot(gs[:,1])

    cmap,norm,ticks = pl.GetDiscreteColormap([-100,0]) # min(vals),0])
    ticks = [float('%1.1f'%(v)) for v in ticks]

    h = ax1.tricontourf(triang,vals,cmap=cmap,norm=norm,levels=ticks,extend='both')

    cb = fig.colorbar(h,cax=cax)#,boundaries=ticks)
    cb.set_label(u'Setzung [mm]',size=12)

    for ke in range(res.nBeams):
        inel = res.beam.inel[ke]
        ax1.plot([-res.coords[2][kn-1] for kn in inel],
                 [res.coords[1][kn-1] for kn in inel],'k',alpha=0.5)

    xlim = [-180,50]
    ax1.plot(xlim,[-4.7,-4.7],'k')
    ax1.plot(xlim,[-9.7,-9.7],'k')
    ax1.plot(xlim,[-31.2,-31.2],'k')

    ax1.set_aspect('equal')
    ax1.set_xlim(xlim)
    ax1.set_ylim([-45,0])

    # cuvette:
    cuv0 = []
    cuv1 = []
    cuv2 = []
    
    shnds = set()
    for ke in range(res.nShells):
        if res.shell.mat[ke] in [13]:
            inel = res.shell.inel[ke]
            for kn in inel:
                shnds.add(kn-1)
    shnds = list(shnds)
    
    for kn in range(res.nNodes):
        if abs(res.coords[1][kn]+0.77)<1e-3 and kn in shnds:
            if abs(res.coords[0][kn])<1e-3:
                cuv0.append((kn,-res.coords[2][kn]))
            elif abs(res.coords[0][kn]+6)<1e-3:
                cuv1.append((kn,-res.coords[2][kn]))
        elif abs(res.coords[1][kn]+1.2)<1e-3:
            if abs(res.coords[0][kn]+2.666666)<1e-3:
                cuv2.append((kn,-res.coords[2][kn]))
    cuv0 = sorted(cuv0,key=lambda v:v[1])
    cuv1 = sorted(cuv1,key=lambda v:v[1])
    cuv2 = sorted(cuv2,key=lambda v:v[1])
    tass0 = [(step.nodal.disp[1][v[0]]-step0.nodal.disp[1][v[0]])*1000 for v in cuv0]
    tass1 = [(step.nodal.disp[1][v[0]]-step0.nodal.disp[1][v[0]])*1000 for v in cuv1]
    tass2 = [(step.nodal.disp[1][v[0]]-step0.nodal.disp[1][v[0]])*1000 for v in cuv2]
    tassmax = min(min(tass0),min(tass1),min(tass2))
    
    col_Uy = []
    diff = []
    
    for kn in range(res.nNodes):
        if abs(res.coords[1][kn]+0.77)<1e-3 and abs(res.coords[0][kn])<1e-3:
            for i in range(6):
                if abs(res.coords[2][kn]-12*i)<1e-3:
                    col_Uy.append((kn,res.coords[2][kn]))
            for i in range(6):
                if abs(res.coords[2][kn]-(66+12*i))<1e-3:
                    col_Uy.append((kn,res.coords[2][kn]))
                
    for i in range(len(col_Uy)-1):
        Uy1 = (step.nodal.disp[1][col_Uy[i+1][0]] - step0.nodal.disp[1][col_Uy[i+1][0]]) * 1000
        Uy0 = (step.nodal.disp[1][col_Uy[i][0]] - step0.nodal.disp[1][col_Uy[i][0]]) * 1000
        diff.append(Uy1-Uy0)
    
    for i in range(len(diff)):
        if i ==(len(diff)-1): 
            ax1.annotate(u'Diff. entre colonnes : %0.1f mm'%(diff[i]),xy=((-col_Uy[i][1]-1),0),#xycoords='axes fraction',
                         ha='right',va='bottom',size=10)
        elif i==5:
            ax1.annotate(u'%0.1f mm'%(diff[i]),xy=((-col_Uy[i][1]-3),0),#xycoords='axes fraction',
                         ha='center',va='bottom',size=10)
        else:
            ax1.annotate(u'%0.1f mm'%(diff[i]),xy=((-col_Uy[i][1]-6),0),#xycoords='axes fraction',
                         ha='center',va='bottom',size=10)
        

    ax0.plot([v[1] for v in cuv0],
             tass0,
             label=u'Bodenplatte Stützenachse')
    ax0.plot([v[1] for v in cuv1],
             [(step.nodal.disp[1][v[0]]-step0.nodal.disp[1][v[0]])*1000 for v in cuv1],
             label=u'Bodenplatte Zwischenachse')
    ax0.plot([v[1] for v in cuv2],
             [(step.nodal.disp[1][v[0]]-step0.nodal.disp[1][v[0]])*1000 for v in cuv2],
             label=u'OK Terrain (435.2 MüM)')

    ax0.set_xlim(xlim)
    ax0.set_ylim([-100,0])
    ax0.set_ylabel('[mm]')
    ax0.annotate(u'Setzung aufgrund Gebäudelasten (T3-T2)\nMax.: %1.1f mm'%(tassmax),xy=(0.5,0.95),xycoords='axes fraction',
                 ha='center',va='top',size=12)
    ax0.legend()
    ax0.grid('both')

    fig.tight_layout()
    fig.savefig(prob+'_tass')
    plt.close(fig)
