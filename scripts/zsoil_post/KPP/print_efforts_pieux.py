# -*- coding: cp1252 -*-
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pickle
from matplotlib.patches import Circle
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools
import time

pathname = '..'
modele = 'M1013_O_2024_v0_5_v24'
pblist = [#modele+'_flach',
##          modele,
##          modele+'_mitGelb',
          modele+'_piles+10']
##          modele+'_piles+10_MEpess']
##          modele+'_piles+10_DN900',
##          modele+'_pilesv1_1',
##          modele+'_pilesv1_2',
##          modele+'_Bopla3m']


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

        tvect = [0,3,6]
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

    step0 = res.steps[tsteps[1]]
    step = res.steps[tsteps[2]]
        
    corner0 = np.array([12.1,-45.97])
    corner1 = np.array([42.37,-7.3])
    corner0 = np.array([70.11,-79.23])
    corner1 = np.array([96.91,-35.43])

    pilecrds = []
    FN = []
    dy_cap = []
    for pile in res.piles:
        beams = []
        fn = []
        for ke0 in pile.beams:
            ke = res.num_beams.index(ke0)
            beams.append(ke)
            inel = res.beam.inel[ke]
            fn.append(step.beam.force[0][ke]-step0.beam.force[0][ke])
        pcrd = (res.coords[0][inel[0]-1],-res.coords[2][inel[0]-1])
        pilecrds.append(pcrd)
        kntop = res.beam.inel[beams[0]][0]
        FN.append(fn)
        # get pile cap settlement:
        dy_cap.append(min(0,step.nodal.disp[1][kntop-1]
                          -step0.nodal.disp[1][kntop-1]))

    # Aufteilung Platte - Pfähle:
    FPiles = sum([f[0] for f in FN])
    sy = 0
    sy0 = 0
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
                if np.mean(x)>corner0[0] and np.mean(x)<corner1[0]:
                    if np.mean(z)>corner0[1] and np.mean(z)<corner1[1]:
                        sy0 += step.cnt.stress[2][ke][kgp]*a*0.25


    shellnodes = list(set([kn-1 for inel in res.shell.inel for kn in inel]))
    bounds = [[min([res.coords[0][kn] for kn in shellnodes]),
               max([res.coords[0][kn] for kn in shellnodes])],
              [min([-res.coords[2][kn] for kn in shellnodes]),
               max([-res.coords[2][kn] for kn in shellnodes])]]
    bounds_glob = bounds
    yref = 397.9
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
        uy.append(min(0,step.nodal.disp[1][kn]
                      -step0.nodal.disp[1][kn]))
        x.append(res.coords[0][kn])
        z.append(-res.coords[2][kn])
    tri = mtri.Triangulation(x,z)
    f = mtri.LinearTriInterpolator(tri,uy)
    xnew = np.arange(bounds_glob[0][0],bounds_glob[0][1],0.5)
    znew = np.arange(bounds_glob[1][0],bounds_glob[1][1],0.5)
    X,Z = np.meshgrid(xnew,znew)
    V = np.zeros([len(X),len(X[0])])
    for kx in range(len(X[0])):
        for kz in range(len(X)):
##            V[kz][kx] = interpolate.griddata((x,z),v,(X[kz][kx],Z[kz][kx]),method='linear')
            V[kz][kx] = f(X[kz][kx],Z[kz][kx])*1000
   
    fig = plt.figure(figsize=(12,14))
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    fig.text(0.01,0.01,prob,size=8)

    ax = fig.add_subplot(111)

    pshells = []
    thick = []
    for ke in range(res.nShells):
        inel = res.shell.inel[ke]
        xy = [np.array([res.coords[0][kn-1],-res.coords[2][kn-1]]) for kn in inel]
        # hatch settlement gap:
        if res.shell.mat[ke]==7:
            pshells.append(Polygon(np.array(xy)))
            if res.shell.EF[ke]==4:
                ax.add_patch(Polygon(np.array(xy),fill=False,hatch='/',lw=0.2))
            thick.append(step.shell.thick[ke])
        # plot wall outlines:
        if res.shell.mat[ke]==8:
            inel = res.shell.inel[ke][:2]
            yy = [res.coords[1][kn-1] for kn in inel]
            if min(yy)<yref+1e-3:
                ax.plot([res.coords[0][kn-1] for kn in inel],
                        [-res.coords[2][kn-1] for kn in inel],'k',
                        lw=1,solid_capstyle='butt')
    polys, cvect = vtktools.get_joined_polygons(pshells, thick)
    pc = PatchCollection(polys,edgecolors='None',
                         cmap=plt.colormaps['Set1'],alpha=0.2)
    pc.set_clim(vmin=min(cvect),vmax=max(cvect))
    pc.set_array(np.array(cvect))
    ax.add_collection(pc)
    
    patches = []
    dscale = 1
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    bbox_props = dict(boxstyle="round,pad=0", fc="w", ec="k", lw=0.4)
    diam = {13:1}
    for k,pile in enumerate(res.piles):
        y = []
        mantel = 0
        L = 0
        for kke,ke0 in enumerate(pile.beams):
            ke = res.num_beams.index(ke0)
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
        mat = res.beam.mat[res.num_beams.index(pile.beams[0])]
        plength = max(y)-min(y)
        patches.append(Circle(pilecrds[k],radius=0.5*diam[mat]/dscale))
        ax.annotate('L%1.0f, %1.0f\n%1.0f'%\
                    (plength,dy_cap[k]*1000,FN[k][0]),
                    xy=pilecrds[k],xytext=(pilecrds[k][0],pilecrds[k][1]),
                    ha='left',va='center',size=8,
                    bbox=bbox_props,rotation=30)
    ax.annotate(u'Länge [m], Setzung Pfahlkopf [mm]\nNormalkraft Pfahlkopf [kN]',
                xy=(55+1,90-12),
                ha='left',va='bottom',size=14,
                bbox=bbox_props,rotation=30)
    Ftop = [FN[k][0] for k in range(len(pilecrds))]
    cmap = plt.colormaps['jet']
    if len(Ftop):
        pc = PatchCollection(patches,edgecolors=('k',),cmap=cmap)
        pc.set_clim(vmin=min(Ftop),vmax=max(Ftop))
        pc.set_array(np.array(Ftop))
        ax.add_collection(pc)
    ax.annotate('Maximale Setzung: %1.0f mm\nLastsumme Bopla 1.5m: '%(min(uy)*1e3)+
                u'%1.1f MN (Platte %1.0f%%, Pfähle %1.0f%%)'%(1e-3*(sy0+sum(Ftop)),
                                                              100*sy0/(sy0+sum(Ftop)),
                                                              100*sum(Ftop)/(sy0+sum(Ftop))),
                xy=(56,23),va='bottom')
    print(sy,sy0)

##    esp = 10.0
##    levels = [int(min(dy_cap)*1000)+k*esp
##              for k in range(int(1000*(max(dy_cap)-min(dy_cap))/esp)+1)]
    CS = ax.contour(X,Z,V,10)
##    CS = ax.contour(X,Z,V,N,linewidths=[(2,) for k in range(N)])
    ax.clabel(CS,inline=0,inline_spacing=2,fontsize=14,colors='k',fmt='%1.0f')
    
    ax.set_xlim(bounds[0])
    ax.set_ylim(bounds[1])
    ax.set_xlim([55,110])
    ax.set_ylim([22,90])
    ax.set_aspect('equal')

    fig.tight_layout()
##    fig.savefig(figname+'_pileres')
    fig.savefig(figname+'_pileres.svg')
    plt.close(fig)

   








