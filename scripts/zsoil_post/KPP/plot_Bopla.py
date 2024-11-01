# -*- coding: cp1252 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pickle
from matplotlib.patches import Circle
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools
from zsoil_tools import postpro_lib as pl
import vtk

pathname = '../..'
modele = 'M1013_O_2024_v0_5_v24'
pblist = [modele+'_piles+10']
t0 = 3
tvect = [3,6]

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

        tsteps = []
        for kt,step in enumerate(res.steps):
            if step.time in tvect and step.sf==0 and step.conv_status==-1:
                tsteps.append(kt)
                
        res.out_steps = tsteps
        res.read_dat()
    ##        res.read_s02()
        res.read_s00()
        res.read_s04()
        res.read_s07()
        pickle.dump(res, open(prob+'.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    step0 = res.steps[tsteps[0]]
    step = res.steps[tsteps[-1]]

    corner0 = np.array([70.11,-79.23])
    corner1 = np.array([96.91,-35.43])

    xlim = [55,110]
    ylim = [22,87]
    y0 = 397.9

    pilecrds = []
    blists = []
    FN = []
    mat = []
    # get pile cap settlement:
    dy_cap = []
    for kp,pile in enumerate(res.piles):
        blist = [res.num_beams.index(ke) for ke in pile.beams]
        blists.append(blist)
        inel = res.beam.inel[blist[0]]
        pcrd = (res.coords[0][inel[0]-1],-res.coords[2][inel[0]-1])
        pilecrds.append(pcrd)
        mat.append(res.beam.mat[blist[0]])

        kntop = res.beam.inel[blist[0]][0]-1
        # get pile cap settlement:
        dy_cap.append(min(0,step.nodal.disp[1][kntop]
                          -step0.nodal.disp[1][kntop]))

        # get normal forces of piles:
        fn = [(step.beam.force[0][ke]-step0.beam.force[0][ke]) for ke in blist]
        FN.append(fn)

    # Aufteilung Platte - Pfähle und Bettungsmoduln:    
    fig = plt.figure(figsize=(7,9))
    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
    fig.text(0.01,0.01,prob,size=8)
    patches = []

    FPiles = sum([f[0] for f in FN])
    sy = 0
    A = 0
    dy_cnt = []
    dp_cnt = []
    p_cnt = []
    dy_cnt_t = []
    dp_cnt_t = []
    for ke in range(res.nContacts):
        if res.cnt.type[ke]==1:
            inel = res.cnt.inel[ke][:4]
            xy = [np.array([res.coords[0][kn-1],
                            -res.coords[2][kn-1]]) for kn in inel]
            if abs(max([res.coords[1][kn-1] for kn in inel])-y0)<0.01:
                poly = Polygon(np.array(xy))
                patches.append(poly)
                dy_cnt.append(0.25*sum([step.nodal.disp[1][kn-1]
                                        -step0.nodal.disp[1][kn-1]for kn in inel]))
                dp_cnt.append(0.25*sum([step.cnt.stress[2][ke][kgp]
                                        -step0.cnt.stress[2][ke][kgp] for kgp in range(4)]))
                p_cnt.append(0.25*sum([step.cnt.stress[2][ke][kgp] for kgp in range(4)]))

                dy_cnt_t.append([0.25*sum([step.nodal.disp[1][kn-1]
                                           -step0.nodal.disp[1][kn-1]for kn in inel]) for kt in tsteps])
                dp_cnt_t.append([0.25*sum([step.cnt.stress[2][ke][kgp]
                                           -step0.cnt.stress[2][ke][kgp] for kgp in range(4)]) for kt in tsteps])

                a = res.compute_area(inel)
                A += a
                # get load distribution only for main slab:
                if np.mean(xy,0)[0]>corner0[0] and np.mean(xy,0)[0]<corner1[0]:
                    if -np.mean(xy,0)[1]>corner0[1] and -np.mean(xy,0)[1]<corner1[1]:
                        for kgp in range(4):
                            sy += (step.cnt.stress[2][ke][kgp])*a*0.25
    dk_cnt = [dp_cnt[k]/dy_cnt[k] for k in range(len(dy_cnt))]
    Gtot = FPiles+sy
    print('Gesamtlast: %1.0f, Pfähle: %1.0f%%, Sohlpressung: %1.0f%%'%
          (Gtot,FPiles/Gtot*100,sy/Gtot*100))

    # Get wall outlines for illustration
    walls = []
    for ke in range(res.nShells):
        if res.shell.mat[ke] in [8]:
            inel = res.shell.inel[ke]
            yy = [res.coords[1][kn-1] for kn in inel]
            if min(yy)<y0+0.01:
                xx = [res.coords[0][kn-1] for kn in inel]
                zz = [-res.coords[2][kn-1] for kn in inel]
                walls.append([[min(xx),max(xx)],[min(zz),max(zz)]])

    ######################################################
    #   Bettung:
    ######################################################
    cmap,norm,ticks = pl.GetDiscreteColormap([min(dk_cnt),max(dk_cnt)])
    cmap,norm,ticks = pl.GetDiscreteColormap([0,5000])
    pc = PatchCollection(patches,edgecolors=('none',),cmap=cmap,norm=norm)
    pc.set_array(np.array(dk_cnt))
    ax = fig.add_subplot(111)
    ax.add_collection(pc)
    cb = fig.colorbar(pc,orientation='horizontal',pad=0.05)
    cb.set_label(u'Bettung verteilt [kN/m3]',size=16)
    cb.ax.tick_params(labelsize=16)

    diamdict = {13:1,14:0.9}

    for wall in walls:
        ax.plot(wall[0],wall[1],'k',lw=1,alpha=0.5,solid_capstyle='butt')

    circles = []
    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=1)
    of = open(prob+'_federkonstanten.csv','w')
    of.write('Pfahlnummer;X [m];Y [m];D [m];L [m];Normalkraft Pfahlkopf [kN];Einsenkung aufgrund Gebäudelast [mm];Federsteifigkeit [MN/m]\n')
    for kp,pile in enumerate(res.piles):
        L = 0
        for ke in blists[kp]:
            inel = res.beam.inel[ke]
            yy = [res.coords[1][kn-1] for kn in inel]
            L += abs(yy[0]-yy[1])
        of.write('%s;%1.2f;%1.2f;%1.1f;%1.0f;%1.0f;%1.3f;%1.3f\n'%
                 (kp+1,pilecrds[kp][0],pilecrds[kp][1],diamdict[mat[kp]],L,FN[kp][0],dy_cap[kp]*1000,FN[kp][0]/dy_cap[kp]/1e3))

        circles.append(Circle(pilecrds[kp],radius=diamdict[mat[kp]]/2.))
        ax.annotate('%1.0f'%(FN[kp][0]/dy_cap[kp]*0.001),
                    xy=pilecrds[kp],xytext=(pilecrds[kp][0]+0.5,pilecrds[kp][1]),
                    ha='left',va='center',size=6)
##                    bbox=bbox_props,rotation=0)
    of.write('\nSumme Normalkraft oberste Pfahlelemente [kN];%1.0f\n'%(FPiles))
    of.write('Summe Sohlpressung [kN];%1.0f\n'%(sy))
    of.write('Lastanteil Pfähle [%%];%1.0f\n'%(FPiles/Gtot*100))
    of.write('Lastanteil Platte [%%];%1.0f\n'%(sy/Gtot*100))
    of.close()

    pc = PatchCollection(circles,edgecolors=('k',),facecolors=('none',))
    ax.add_collection(pc)
    ax.annotate('Bettung verteilt in [kN/m3]\nPfahlsteifigkeiten in [MN/m]',
                xy=(57,65),size=10,ha='left',va='bottom',rotation=90)

    ax.set_aspect('equal')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    fig.tight_layout()
    fig.savefig(prob+'_Bettung')
##    fig.savefig('svg/'+prob+'_Bettung.svg')
   
    ######################################################
    #  Sohlpressung:    
    ######################################################
    fig = plt.figure(figsize=(7,9))
    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
    fig.text(0.01,0.01,prob,size=8)
    patches3 = []
    cvect = []
    cscale = []
    minmax = [[1e10,-1e10],[1e10,-1e10]]


##    # check self-weight of BP:
##    weightBP = 0
##    for ke in range(res.nShells):
##        inel = res.shell.inel[ke]
##        x = [res.coords[0][kn-1] for kn in inel]
##        z = [res.coords[2][kn-1] for kn in inel]
##        a = (max(x)-min(x))*(max(z)-min(z))
##        th = step.shell.thick[ke]
##        weightBP += a*th*25
        
##    cmap,norm,ticks = pl.GetDiscreteColormap([-82.71470069885254,0])
    cmap,norm,ticks = pl.GetDiscreteColormap([min(p_cnt),0])
    if 'flach' in prob:
        cmap,norm,ticks = pl.GetDiscreteColormap([-800,0])
    else:
        cmap,norm,ticks = pl.GetDiscreteColormap([-400,0])
    pc = PatchCollection(patches,edgecolors=('none',),cmap=cmap,norm=norm)
    pc.set_array(np.array(p_cnt))
    ax = fig.add_subplot(111)
    ax.add_collection(pc)
    cb = fig.colorbar(pc,orientation='horizontal',pad=0.05)
    cb.set_label(u'Sohlpressung [kPa]',size=16)
    cb.ax.tick_params(labelsize=16)

    for wall in walls:
        ax.plot(wall[0],wall[1],'k',lw=1,alpha=0.5,solid_capstyle='butt')
    
    for kp in range(res.nPiles):
        ax.add_patch(Circle(pilecrds[kp],radius=diamdict[mat[kp]]/2.,ec='k',fc='None'))
        ax.annotate('%1.0f'%(FN[kp][0]),
                    xy=pilecrds[kp],xytext=(pilecrds[kp][0]+0.5,pilecrds[kp][1]),
                    ha='left',va='center',size=6)

    minp = min(p_cnt)
    min_xy = np.mean(patches[p_cnt.index(minp)].xy,0)
    ax.annotate('%1.0f kPa'%(minp),xy=min_xy,size=14)
    ax.plot(min_xy[0],min_xy[1],'k.')

    ax.annotate(u'Gesamtlast: %1.1f MN\n(Platte: %1.0f%%, Pfähle: %1.0f%%)\nKräfte Pfahlkopf in [kN]'%
                (Gtot*1e-3,sy/Gtot*100,FPiles/Gtot*100),
                xy=(123,115),size=12,ha='right',va='top')

    ax.set_aspect('equal')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    fig.tight_layout()
    fig.savefig(prob+'_Sohlpressung')
##    fig.savefig('svg/'+prob+'_Sohlpressung.svg')
    plt.close()
   

if False:
    # Last-Setzungsverläufe:
    fig = plt.figure(figsize=(11,10))
    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
    fig.text(0.01,0.01,prob,size=8)
        
    ax = fig.add_subplot(111)
    scy=50
    scp=0.002
    for kp in range(1,len(dy_cnt_t),1):
        ax.plot([v*scy+np.mean(patches[kp].xy,0)[0] for v in dy_cnt_t[kp]],
                [v*scp+np.mean(patches[kp].xy,0)[1] for v in dp_cnt_t[kp]],
                'k-',lw=0.5)
    for kp in range(res.nPiles):
        ax.add_patch(Circle(pilecrds[kp],radius=diamdict[mat[kp]]/5.,ec='r',fc='None'))

    ax.set_aspect('equal')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    fig.tight_layout()
    fig.savefig(prob+'_Last-Setzung')




