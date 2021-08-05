# -*- coding: cp1252 -*-
import numpy as np
import math
import matplotlib.pyplot as plt
import vtk
import re

from zsoil_tools import vtktools

def plot_disp_paroi(prob,coupes,etapes,ketapes,niveaux,res,figsize=(7,10),figpath='.'):
    shnds = set()
    shele = set()
    for ke in range(res.nShells):
        if res.shell.mat[ke] in [1,9,30] and res.shell.EF[ke]==1:
            inel = res.shell.inel[ke]
            xm = np.mean([res.coords[0][kn-1] for kn in inel])
            ym = np.mean([res.coords[1][kn-1] for kn in inel])
            zm = np.mean([res.coords[2][kn-1] for kn in inel])
            if ym>691.3:
                if (xm-(26702.03607-26600))**2+(zm-(11800-11619.95921))**2<33.65854**2:
                    if xm>100.414:
                        for kn in inel:
                            shnds.add(kn-1)
                        shele.add(ke)
                else:
                    for kn in inel:
                        shnds.add(kn-1)
                    shele.add(ke)
    shnds = list(shnds)
    shele = list(shele)

    p0 = np.array([101,138.5])
    nlist = []
    anglist = []
    for kn in range(res.nNodes):
        crd = np.array([res.coords[0][kn],res.coords[1][kn],res.coords[2][kn]])
        if kn in shnds:
            ang = np.angle((crd[0]-p0[0]+(crd[2]-p0[1])*1j)*(0.2+1j))
            if ang in anglist:
                nlist[anglist.index(ang)][0].append((kn,crd))
            else:
                nlist.append(([(kn,crd)],ang))
                anglist.append(ang)

    nlist.sort(key=lambda v:v[1])
    anglist = [v[1] for v in nlist]
    nlist = [v[0] for v in nlist]
    xn = [v[0][1][0] for v in nlist]
    zn = [-v[0][1][2] for v in nlist]

    figs = []
    for kniv,niveau in enumerate(niveaux):    
        fig = plt.figure(figsize=(12,8))
        fig.text(0.01,0.01,prob,size=8)
        ax = plt.subplot(111)
        
        for ke in range(res.nTrusses):
            inel = res.truss.inel[ke]
            xx = [res.coords[0][kn-1] for kn in inel]
            yy = [res.coords[1][kn-1] for kn in inel]
            zz = [-res.coords[2][kn-1] for kn in inel]
            if res.truss.EF[ke] in [20]:
                ax.plot(xx,zz,'k',alpha=0.5)
        
    ##    for ke in range(res.nNails):
    ##        beams = res.nails[ke].beams
    ##        inel = [res.beam.inel[res.num_beams.index(beams[0])][0],
    ##                res.beam.inel[res.num_beams.index(beams[-1])][1]]
    ##        xx = [res.coords[0][kn-1] for kn in inel]
    ##        yy = [res.coords[1][kn-1] for kn in inel]
    ##        zz = [-res.coords[2][kn-1] for kn in inel]
    ####        if res.beam.EF[ke] in [18]:
    ##        ax.plot(xx,zz,'k',alpha=0.5)
        
    ##    for ke in range(res.nAnchors):
    ##        trusses = res.anchors[ke].trusses
    ##        inel = [res.truss.inel[res.num_trusses.index(trusses[0])][0],
    ##                res.truss.inel[res.num_trusses.index(trusses[-1])][1]]
    ##        xx = [res.coords[0][kn-1] for kn in inel]
    ##        yy = [res.coords[1][kn-1] for kn in inel]
    ##        zz = [-res.coords[2][kn-1] for kn in inel]
    ####        if res.beam.EF[ke] in [18]:
    ##        ax.plot(xx,zz,'k',alpha=0.5)

        tvect = [etapes.tvect[kk] for kk in ketapes]
        enames = [etapes.names[kk] for kk in ketapes]

        scale = 5.e2
        for kkt,kt in enumerate(res.out_steps):
            step = res.steps[kt]
            if step.time in tvect:
                ind = tvect.index(step.time)
                dhx = []
                dhz = []
                dhmaxlist = []
                for nn in nlist:
                    dhvect = [(step.nodal.disp[0][v[0]]**2+step.nodal.disp[2][v[0]]**2)**0.5 for v in nn]
##                    # déplacement maximal de la paroi:
##                    dhmax = max(dhvect)
##                    # déplacement en tête de la paroi:
##                    yvect = [res.coords[1][kn] for kn in [v[0] for v in nn]]
##                    ymaxind = yvect.index(max(yvect))
##                    dhmax = dhvect[ymaxind]
                    # déplacement le plus proche du niveau x:
                    dyvect = [abs(res.coords[1][kn]-niveau) for kn in [v[0] for v in nn]]
                    yind = dyvect.index(min(dyvect))
                    dhmax = dhvect[yind]
                    dhmaxlist.append(dhmax)
                    knmax = nn[dhvect.index(dhmax)][0]
                    dhx.append(step.nodal.disp[0][knmax])
                    dhz.append(-step.nodal.disp[2][knmax])
                xmaxind = xn.index(max(xn))
                h = ax.plot([xn[kk]+scale*dhx[kk] for kk in range(xmaxind+1)],
                            [zn[kk]+scale*dhz[kk] for kk in range(xmaxind+1)])
                ax.plot([xn[kk]+scale*dhx[kk] for kk in range(xmaxind+1,len(nlist))],
                        [zn[kk]+scale*dhz[kk] for kk in range(xmaxind+1,len(nlist))],
                        color=h[0].get_color(),
                        label='T=%1.0f: %s'%(tvect[ind],enames[ind]))
            if step.time==tvect[-1]:
                parois = [(0.75,1.37), # paroi nord
                          (-1.79,0.64),  # paroi Ouest 1
                          (-2.39,-1.79),   # paroi Ouest 2
                          (1.94,3.02)] # paroi Sud
                for paroi in parois:
                    try:
                        angind0 = next(x[0] for x in enumerate(anglist) if x[1]>paroi[0])
                    except:
                        angind0 = len(anglist)-1
                    try:
                        angind1 = next(x[0] for x in enumerate(anglist) if x[1]>paroi[1])
                    except:
                        angind1 = len(anglist)-1

                    dhmax = max(dhmaxlist[angind0:angind1+1])
                    dhmaxind = dhmaxlist.index(dhmax)
                    xx = [xn[dhmaxind],xn[dhmaxind]+scale*dhx[dhmaxind]]
                    zz = [zn[dhmaxind],zn[dhmaxind]+scale*dhz[dhmaxind]]
                    ax.annotate('',xy=(xx[1],zz[1]),
                                xytext=(xx[0],zz[0]),
                                arrowprops=dict(arrowstyle='-|>'))
                    ax.annotate('$u_h$=%1.0f mm'%(dhmax*1e3),xy=(np.mean(xx),np.mean(zz)),
                                bbox=dict(facecolor='w',edgecolor='k',boxstyle='round,pad=0.2'))
                
    ##        ax.grid('on')
        ax.set_aspect('equal')
        ax.axis('off')
        ax.legend()
        
        for ke in range(res.nShells):
            if res.shell.mat[ke] in [1,9,25,30] and ke in shele:
                inel = res.shell.inel[ke]
                xx = [res.coords[0][kn-1] for kn in inel]
                yy = [res.coords[1][kn-1] for kn in inel]
                zz = [-res.coords[2][kn-1] for kn in inel]
                if abs(min(yy)-691.3)<1e-3:
                    ax.plot(xx,zz,'k')


        for k in range(len(coupes.origins)):
            a0 = -10
            a1 = 10
            loc_syst = np.array([np.cross(coupes.normals[k],(0,1,0)),(0,1,0)])
            ax.plot([coupes.origins[k][0]+a0*loc_syst[0][0],coupes.origins[k][0]+a1*loc_syst[0][0]],
                    [-coupes.origins[k][2]-a0*loc_syst[0][2],-coupes.origins[k][2]-a1*loc_syst[0][2]],
                    'k--',alpha=0.5)
            ax.annotate(coupes.names[k],xy=(coupes.origins[k][0]+a0*loc_syst[0][0],-coupes.origins[k][2]-a0*loc_syst[0][2]))

        if niveau>800:
            ax.annotate(u'Déplacement horizontal en tête de la paroi',xy=(0.1,0.1),xycoords='axes fraction',size=14)
        else:
            ax.annotate(u'Déplacement horizontal à la cote %1.1f msm'%(niveau),xy=(0.1,0.1),xycoords='axes fraction',size=14)

        fig.tight_layout()
        fig.savefig(figpath+'/'+prob+'_%1.0fmsm_disp_paroi'%(niveau))
        plt.close(fig)
        figs.append(figpath+'/'+prob+'_%1.0fmsm_disp_paroi'%(niveau))

    return figs
