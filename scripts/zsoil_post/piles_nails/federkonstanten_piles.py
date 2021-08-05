# -*- coding: cp1252 -*-
import numpy,os,math,cmath
from numpy import linalg as la
from numpy import fft
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pickle
from matplotlib.patches import Circle
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy import interpolate
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from zsoil_tools import zsoil_results as zr
from zsoil_tools import postpro_lib as pl

pathname = '..'
pblist = ['M1230_v1_1']


for kf,prob in enumerate(pblist):
    figname = prob
    try:
        res = pickle.load(open(prob+'_bettung.p', "rb" ))
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
        pickle.dump(res, open(prob+'_bettung.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    plabs = ['C%i'%(k+1) for k in range(res.nPiles)]

    Beams = []
    py = []
    pilecrds = []
    FN = []
    mat = []
    for pile in res.piles:
        beams = []
        py0 = []
        for ke0 in pile.beams:
            ke = res.num_beams.index(ke0)
            beams.append(ke)
            inel = res.beam.inel[ke]
            py0.append(0.5*(res.coords[1][inel[0]-1]+res.coords[1][inel[1]-1]))
        py.append(py0)
        pcrd = (res.coords[0][inel[0]-1],-res.coords[2][inel[0]-1])
        pilecrds.append(pcrd)
        Beams.append(beams)
        ke = res.num_beams.index(pile.beams[0])
        mat.append(res.beam.mat[ke])

    kntop = []  # 1st node from top of pile
    for kp,pile in enumerate(res.piles):
        kntop.append(res.beam.inel[Beams[kp][0]][0])
    # get pile cap settlement:
    step0 = res.steps[tsteps[0]]
    step = res.steps[tsteps[1]]
    dy_cap = []
    for kp in range(res.nPiles):
        dy_cap.append(min(0,step.nodal.disp[1][kntop[kp]-1]
                          -step0.nodal.disp[1][kntop[kp]-1]))
    

    # get normal forces of piles:
    for el in Beams:
        fn = []
        for ke in el:
            fn.append(step.beam.force[0][ke]-step0.beam.force[0][ke])
        FN.append(fn)

    # Aufteilung Platte - Pfähle und Bettungsmoduln:    
    fig = plt.figure(figsize=(16,16))
    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
    fig.text(0.01,0.01,prob,size=8)
    fig.text(0.26,0.02,'Bettung verteilt [kN/m3] und punktuell [MN/m]',size=18,bbox=bbox_props)
    patches = []
    cvect = []
    minmax = [[1e10,-1e10],[1e10,-1e10]]

    FPiles = sum([f[0] for f in FN])
    sy = 0
    A = 0
    for ke in range(res.nContacts):
        if res.cnt.type[ke]==1:
            xy = numpy.zeros(4*2)
            xy.shape = (4,2)
            inel = res.cnt.inel[ke][:4]
            for kk,kn in enumerate(inel):
                xy[kk][0] = res.coords[0][kn-1]
                xy[kk][1] = -res.coords[2][kn-1]
                minmax[0][0] = min(minmax[0][0],xy[kk][0])
                minmax[0][1] = max(minmax[0][1],xy[kk][0])
                minmax[1][0] = min(minmax[1][0],xy[kk][1])
                minmax[1][1] = max(minmax[1][1],xy[kk][1])
            poly = Polygon(numpy.array(xy))
            patches.append(poly)
            dy_cnt = 0.25*sum([step.nodal.disp[1][kn-1]
                               -step0.nodal.disp[1][kn-1]for kn in inel])
            p_cnt = 0.25*sum([step.cnt.stress[2][ke][kgp]
                              -step0.cnt.stress[2][ke][kgp] for kgp in range(4)])
            cvect.append(p_cnt/dy_cnt)
            
            x = [res.coords[0][kn-1] for kn in inel]
            z = [res.coords[2][kn-1] for kn in inel]
            a = (max(x)-min(x))*(max(z)-min(z))
            A += a
            for kgp in range(4):
                sy += (step.cnt.stress[2][ke][kgp]
                       -step0.cnt.stress[2][ke][kgp])*a*0.25
    Gtot = FPiles+sy
    print('Gesamtlast: %1.0f, Pfähle: %1.0f%%, Sohlpressung: %1.0f%%'%
          (Gtot,FPiles/Gtot*100,sy/Gtot*100))

    cmap,norm,ticks = pl.GetDiscreteColormap([-1000,10000])
    pc = PatchCollection(patches,edgecolors=('none',),cmap=cmap,norm=norm)
    pc.set_array(numpy.array(cvect))
    ax = fig.add_subplot(111)
    ax.add_collection(pc)
    ax.set_xlim(minmax[0])
    ax.set_ylim(minmax[1])
    cb = fig.colorbar(pc)
    cb.set_label(u'Bettung verteilt [kN/m3]',size=16)
    cb.ax.tick_params(labelsize=16)

    diamdict = {11:0.6,12:0.9}

    patches = []
    of = open(prob+'_federkonstanten.csv','w')
    of.write('Pfahlnummer;X [m];Y [m];D [m];Normalkraft Pfahlkopf [kN];Einsenkung aufgrund Gebäudelast [mm]\n')
    for kp in range(res.nPiles):
        of.write('%s;%1.2f;%1.2f;%1.1f;%1.0f;%1.3f\n'%
                 (plabs[kp],pilecrds[kp][0],pilecrds[kp][1],diamdict[mat[kp]],FN[kp][0],dy_cap[kp]*1000))
        
        patches.append(Circle(pilecrds[kp],radius=diamdict[mat[kp]]/2.))
        ax.annotate('%1.0f'%(FN[kp][0]/dy_cap[kp]*0.001),
                    xy=pilecrds[kp],xytext=(pilecrds[kp][0],pilecrds[kp][1]),
                    ha='left',va='center',size=9,
                    bbox=bbox_props,rotation=0)
    of.write('\nSumme Normalkraft oberste Pfahlelemente [kN];%1.0f\n'%(FPiles))
    of.write('Summe Sohlpressung [kN];%1.0f\n'%(sy))
    of.write('Lastanteil Pfähle [%%];%1.0f\n'%(FPiles/Gtot*100))
    of.write('Lastanteil Platte [%%];%1.0f\n'%(sy/Gtot*100))
    of.close()

    pc = PatchCollection(patches,edgecolors=('k',),facecolors=('none',))
    ax.add_collection(pc)

    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(prob+'_Bettung')
   
    # Sohlpressung:    
    fig = plt.figure(figsize=(16,16))
    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
    fig.text(0.01,0.01,prob,size=8)
    patches = []
    cvect = []
    minmax = [[1e10,-1e10],[1e10,-1e10]]

    FPiles = sum([f[0] for f in FN])
    sy = 0
    A = 0
    for ke in range(res.nContacts):
        if res.cnt.type[ke]==1:
            xy = numpy.zeros(4*2)
            xy.shape = (4,2)
            inel = res.cnt.inel[ke][:4]
            for kk,kn in enumerate(inel):
                xy[kk][0] = res.coords[0][kn-1]
                xy[kk][1] = -res.coords[2][kn-1]
                minmax[0][0] = min(minmax[0][0],xy[kk][0])
                minmax[0][1] = max(minmax[0][1],xy[kk][0])
                minmax[1][0] = min(minmax[1][0],xy[kk][1])
                minmax[1][1] = max(minmax[1][1],xy[kk][1])
            poly = Polygon(numpy.array(xy))
            patches.append(poly)
            p_cnt = 0.25*sum([step.cnt.stress[2][ke][kgp]
                              -step0.cnt.stress[2][ke][kgp] for kgp in range(4)])
            cvect.append(p_cnt)
            
            x = [res.coords[0][kn-1] for kn in inel]
            z = [res.coords[2][kn-1] for kn in inel]
            a = (max(x)-min(x))*(max(z)-min(z))
            A += a
            for kgp in range(4):
                sy += (step.cnt.stress[2][ke][kgp]
                       -step0.cnt.stress[2][ke][kgp])*a*0.25
    Gtot = FPiles+sy

##    # check self-weight of BP:
##    weightBP = 0
##    for ke in range(res.nShells):
##        inel = res.shell.inel[ke]
##        x = [res.coords[0][kn-1] for kn in inel]
##        z = [res.coords[2][kn-1] for kn in inel]
##        a = (max(x)-min(x))*(max(z)-min(z))
##        th = step.shell.thick[ke]
##        weightBP += a*th*25
        
    cmap,norm,ticks = pl.GetDiscreteColormap([-60,0])
    pc = PatchCollection(patches,edgecolors=('none',),cmap=cmap,norm=norm)
    pc.set_array(numpy.array(cvect))
    ax = fig.add_subplot(111)
    ax.add_collection(pc)
    ax.set_xlim(minmax[0])
    ax.set_ylim(minmax[1])
    cb = fig.colorbar(pc)
    cb.set_label(u'Sohlpressung [kPa]',size=16)
    cb.ax.tick_params(labelsize=16)

    fig.text(0.5,0.04,u'Gesamtlast: %1.1f MN (Platte: %1.0f%%, Pfähle: %1.0f%%)'%
             (Gtot*1e-3,sy/Gtot*100,FPiles/Gtot*100),
             size=16,horizontalalignment='center')

    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(prob+'_Sohlpressung')
    plt.close()
   








