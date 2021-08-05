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
from openpyxl import load_workbook
from scipy import interpolate
import six
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from zsoil_tools import zsoil_results as zr
import time

pathname = '..'
##pathname = '//192.168.1.51/Mandats sur F RAID0/M1013_Hardturm'
pblist = ['M1013_O2021_v1_pilesV1017_cent1_L1Ost',
          'M1013_O2021_v1_pilesV1017_cent1_L2Ost',
          'M1013_O2021_v1_piles_red1_L1Ost',
          'M1013_O2021_v1_noPiles_L1Ost']
pblist = ['M1013_O_v8_pilesV1017_cent2_l171110']


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

        tvect = [0,4,5]
        tsteps = []
        for kt,step in enumerate(res.steps):
            if step.time in tvect and step.sf==0 and step.conv_status==-1:
                tsteps.append(kt)
##        tsteps = [5,6]
                
        res.out_steps = tsteps
        res.read_dat()
        res.read_s02()
        res.read_s00()
        res.read_s04()
        res.read_s07()
        pickle.dump(res, open(prob+'.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    plabs = ['C%i'%(k) for k in range(1,84)]

    Beams = []
    py = []
    pilecrds = []
    FN = []
    for pile in res.piles:
        beams = []
        py0 = []
        for ke0 in pile.beams:
            ke = res.num_beams.index(ke0)
            beams.append(ke)
            inel = res.beam.inel[ke]
            py0.append(0.5*(res.coords[1][inel[0]-1]
                            +res.coords[1][inel[1]-1]))
        py.append(py0)
        pcrd = (res.coords[0][inel[0]-1],-res.coords[2][inel[0]-1])
        pilecrds.append(pcrd)
        Beams.append(beams)


    kntop = []  # 1st node from top of pile
    for kp,pile in enumerate(res.piles):
        kntop.append(res.beam.inel[Beams[kp][0]][0])
    # get pile cap settlement:
    step0 = res.steps[tsteps[1]]
    step = res.steps[tsteps[2]]
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

    # Aufteilung Platte - Pfähle:
    FPiles = sum([f[0] for f in FN])
    sy = 0
    A = 0
    for ke in range(res.nContacts):
        if res.cnt.type[ke]==1:
            inel = res.cnt.inel[ke][:4]
            x = [res.coords[0][kn-1] for kn in inel]
            z = [res.coords[2][kn-1] for kn in inel]
            a = (max(x)-min(x))*(max(z)-min(z))
            A += a
            for kgp in range(4):
                sy += step.cnt.stress[2][ke][kgp]*a*0.25


##    bounds = [[min([crd[0] for crd in pilecrds])-2,
##               max([crd[0] for crd in pilecrds])+2],
##              [min([crd[1] for crd in pilecrds])-2,
##               max([crd[1] for crd in pilecrds])+2]]
##    bounds_glob = [[-269, -220], [269, 324]]
##    bounds = [[-269, -220], [269, 324]]
    shellnodes = list(set([kn-1 for inel in res.shell.inel for kn in inel]))
    bounds = [[min([res.coords[0][kn] for kn in shellnodes]),
               max([res.coords[0][kn] for kn in shellnodes])],
              [min([-res.coords[2][kn] for kn in shellnodes]),
               max([-res.coords[2][kn] for kn in shellnodes])]]
    bounds_glob = bounds
    yref = 397
    # contour of setzungen:
##    get nodes of Bodenplatte
    nlist = set()
    for ke in range(res.nVolumics):
        for kn in res.vol.inel[ke]:
            y = res.coords[1][kn-1]
            if abs(y-yref)<1e-6:
                nlist.add(kn-1)
##    for ke in range(res.nShells):
##        if res.shell.mat[ke]==7:
##            for kn in res.shell.inel[ke]:
##                nlist.add(kn-1)
    nlist = list(nlist)
    x = []
    z = []
    uy = []
    for kn in nlist:
        uy.append(min(0,res.steps[tsteps[-1]].nodal.disp[1][kn]
                      -res.steps[tsteps[0]].nodal.disp[1][kn]))
        x.append(res.coords[0][kn])
        z.append(-res.coords[2][kn])
    tri = mtri.Triangulation(x,z)
    f = mtri.LinearTriInterpolator(tri,uy)
##    f = interpolate.Rbf(x,z,v,function='linear')
##    f = interpolate.Rbf(x,z,v,function='linear')
    xnew = numpy.arange(bounds_glob[0][0],bounds_glob[0][1],0.5)
    znew = numpy.arange(bounds_glob[1][0],bounds_glob[1][1],0.5)
    X,Z = numpy.meshgrid(xnew,znew)
    V = numpy.zeros([len(X),len(X[0])])
    for kx in range(len(X[0])):
        for kz in range(len(X)):
##            V[kz][kx] = interpolate.griddata((x,z),v,(X[kz][kx],Z[kz][kx]),method='linear')
            V[kz][kx] = f(X[kz][kx],Z[kz][kx])*1000
   
    fig = plt.figure(figsize=(14,18))
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    fig.text(0.01,0.01,prob,size=8)

    ax = fig.add_subplot(111)

    pshells = []
    thick = []
    for ke in range(res.nShells):
        xy = numpy.zeros(4*2)
        xy.shape = (4,2)
        for kk,kn in enumerate(res.shell.inel[ke]):
            xy[kk] = [res.coords[0][kn-1],-res.coords[2][kn-1]]
        pshells.append(Polygon(numpy.array(xy)))
        if res.shell.EF[ke]==4:
            ax.add_patch(Polygon(numpy.array(xy),fill=False,hatch='/',lw=0.2))
        thick.append(step.shell.thick[ke])
    pc = PatchCollection(pshells,edgecolors=('None',),
                         cmap=plt.cm.get_cmap('Set1'),alpha=0.2)
    pc.set_clim(vmin=min(thick),vmax=max(thick))
    pc.set_array(numpy.array(thick))
    ax.add_collection(pc)

    for ke in range(res.nBeams):
        if res.beam.mat[ke]==8:
            inel = res.beam.inel[ke]
            ax.plot([res.coords[0][kn-1] for kn in inel],
                    [-res.coords[2][kn-1] for kn in inel],'k')
    
    patches = []
    dscale = 1
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    bbox_props = dict(boxstyle="round,pad=0", fc="w", ec="k", lw=0.4)
    diam = {13:0.4,14:0.5}
    for k in range(res.nPiles):
        el = Beams[k]
        y = []
        mantel = 0
        L = 0
        for kke,ke in enumerate(el):
            y.append(res.coords[1][res.beam.inel[ke][0]-1])
            y.append(res.coords[1][res.beam.inel[ke][1]-1])
            dl = abs(res.coords[1][res.beam.inel[ke][0]-1]
                     -res.coords[1][res.beam.inel[ke][1]-1])
            L += dl
            mantel += 0.5*sum([v*dl for v in step.cnt.str_level[res.num_contacts.index(res.piles[k].cnt1D[kke])]])
##        print('Pile %i: Mobilized sleeve friction %1.1f%%'%(k+1,mantel/L*100))
##        print('Pile %i: Mobilized tip resistance %1.1f%%'%\
##              (k+1,100*step.cnt.str_level\
##               [res.num_contacts.index(res.piles[k].cnt0D)][0]))
        mat = res.beam.mat[el[0]]
        plength = max(y)-min(y)
        patches.append(Circle(pilecrds[k],radius=diam[mat]/dscale))
        ax.annotate(plabs[k]+', L%1.0fm\n%1.0fmm\n%1.0fkN, %1.0fkN'%\
                    (plength,dy_cap[k]*1000,FN[k][0],FN[k][-1]),
                    xy=pilecrds[k],xytext=(pilecrds[k][0],pilecrds[k][1]),
                    ha='left',va='center',size=8,
                    bbox=bbox_props,rotation=45)
    ax.annotate(u'Pfahlnummer, Länge\nSetzung Pfahlkopf\nNormalkraft Pfahlkopf, -spitze',
                xy=(bounds[0][0]+1,bounds[1][1]-11),
                ha='left',va='bottom',size=16,
                bbox=bbox_props,rotation=45)
    Ftop = [FN[k][0] for k in range(len(pilecrds))]
    cmap = plt.cm.get_cmap('jet')
    if len(Ftop):
        pc = PatchCollection(patches,edgecolors=('k',),cmap=cmap)
        pc.set_clim(vmin=min(Ftop),vmax=max(Ftop))
        pc.set_array(numpy.array(Ftop))
        ax.add_collection(pc)

##    esp = 10.0
##    levels = [int(min(dy_cap)*1000)+k*esp
##              for k in range(int(1000*(max(dy_cap)-min(dy_cap))/esp)+1)]
    CS = ax.contour(X,Z,V,10)
##    CS = ax.contour(X,Z,V,N,linewidths=[(2,) for k in range(N)])
    ax.clabel(CS,inline=0,inline_spacing=2,fontsize=14,colors='k',fmt='%1.0f')

##    ax.set_xticks([53.12+k*6.7 for k in range(7)])
##    ax.set_xticklabels(['A','B','C','D','E','F','G'])
##    ax.set_yticks([117.73-k*7.3 for k in range(7)])
##    ax.set_yticklabels([str(6+k*5) for k in range(7)])    
    
    ax.set_xlim(bounds[0])
    ax.set_ylim(bounds[1])
    ax.set_aspect('equal')
    
##    for x in ax.get_xticks():
##        ax.plot([x,x],bounds[1],'k:')
##    for y in ax.get_yticks():
##        ax.plot(bounds[0],[y,y],'k:')

    fig.tight_layout()
    fig.savefig(figname+'_pileres')
    fig.savefig(figname+'_pileres.svg')
    plt.close(fig)

   








