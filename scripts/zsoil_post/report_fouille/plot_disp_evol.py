# -*- coding: cp1252 -*-
import numpy as np
import math
import matplotlib.pyplot as plt
import pickle


##from zsoil_tools import zsoil_results as zr
##
##tlabvect = [3,5,7,8,9,10,13]
##
##etapes = [u'T=3: exc. 705 msm',
##          u'T=5: exc. 703 msm',
##          u'T=7: exc. 698.5 msm',
##          u'T=8: exc. 696.25 msm',
##          u'T=9: exc. 694 msm',
##          u'T=10: exc. 691.3 msm, clous sud',
##          u'T=13: f.f. atteint']

def plot_disp_evol(prob,etapes,res,figsize=(10,10),figpath='.'):
    ngroups = [u'Paroi devant BIO4',
               u'Paroi Ouest étayée',
               u'Paroi Ouest ancrée 1',
               u'Paroi Ouest ancrée 2',
               u'Paroi Sud',
               u'M2',
               u'Semelle BIO4']
    nlistdefs = [([1,30],[(143.5,163.3),(0,1000),(121.4,131.5)]),
                 ([1,30],[(103.0,127.1),(0,1000),(113.2,130.6)]),
                 ([1,30],[(93.3,100.5),(0,1000),(133.8,139.2)]),
                 ([1,30],[(80.4,92.1),(0,1000),(140.6,154.3)]),
                 ([1,30],[(107.6,146.9),(0,1000),(150.0,169.7)]),
                 ([12],[(80.3,131.8),(0,1000),(96.4,133.1)]),
                 ([24],[(144.2,157.6),(0,1000),(120.0,126.8)])]
##                 ([24],[(144.2,164.6),(0,1000),(120.0,129.4)])]
    nlists = [[] for k in ngroups]

    shnds = [set(),set(),set()]
    for ke in range(res.nShells):
        if res.shell.mat[ke] in [1,30] and res.shell.EF[ke]==1:
            inel = res.shell.inel[ke]
            for kn in inel:
                shnds[0].add(kn-1)
        elif res.shell.mat[ke]==12:
            inel = res.shell.inel[ke]
            for kn in inel:
                shnds[1].add(kn-1)
        elif res.shell.mat[ke]==24:
            inel = res.shell.inel[ke]
            for kn in inel:
                shnds[2].add(kn-1)

    for kn in range(res.nNodes):
        crd = np.array([res.coords[0][kn],res.coords[1][kn],res.coords[2][kn]])
        if kn in shnds[0]:
            for kl in [0,1,2,3,4]:
                ndef = nlistdefs[kl]
                if crd[0]>ndef[1][0][0] and crd[0]<ndef[1][0][1]:
                    if crd[1]>ndef[1][1][0] and crd[1]<ndef[1][1][1]:
                        if crd[2]>ndef[1][2][0] and crd[2]<ndef[1][2][1]:
                            nlists[kl].append(kn)
        elif kn in shnds[1]:
            for kl in [5]:
                ndef = nlistdefs[kl]
                if crd[0]>ndef[1][0][0] and crd[0]<ndef[1][0][1]:
                    if crd[1]>ndef[1][1][0] and crd[1]<ndef[1][1][1]:
                        if crd[2]>ndef[1][2][0] and crd[2]<ndef[1][2][1]:
                            nlists[kl].append(kn)
        elif kn in shnds[2]:
            for kl in [6]:
                ndef = nlistdefs[kl]
                if crd[0]>ndef[1][0][0] and crd[0]<ndef[1][0][1]:
                    if crd[1]>ndef[1][1][0] and crd[1]<ndef[1][1][1]:
                        if crd[2]>ndef[1][2][0] and crd[2]<ndef[1][2][1]:
                            nlists[kl].append(kn)

    Vmin = [[] for kl in ngroups]
    Hmax = [[] for kl in ngroups]
    tx = []
    for kkt,kt in enumerate(res.out_steps):
        step = res.steps[kt]
        tx.append(step.time)
        for kl in range(len(ngroups)):
            hmax = -1e10
            vmin = 1e10
            for kn in nlists[kl]:
                hmax = max(hmax,(step.nodal.disp[0][kn]**2+step.nodal.disp[2][kn]**2)**0.5)
                vmin = min(vmin,step.nodal.disp[1][kn])
            Hmax[kl].append(hmax*1e3)
            Vmin[kl].append(vmin*1e3)
    
    fig = plt.figure(figsize=figsize)
    fig.text(0.01,0.01,prob,size=8)
    AX = [plt.subplot(2,1,1),
          plt.subplot(2,1,2)]

    for kl in range(len(ngroups)):
        AX[0].plot(tx,Hmax[kl],label=ngroups[kl])
        if not 'Paroi' in ngroups[kl]:
            AX[1].plot(tx,Vmin[kl],label=ngroups[kl])

    for ax in AX:
##        ax.set_xlabel(u'Etape [-]',size=14)
        ax.legend()
        ax.set_xticks(etapes.tvect)
        ax.set_xticklabels(['T=%1.0f: %s'%(etapes.tvect[kk],etapes.names[kk]) for kk in range(len(etapes.tvect))])
        plt.setp(ax.xaxis.get_majorticklabels(),rotation=30,ha='right')
        ax.set_xlim([0,max(etapes.tvect)])
        ax.grid('on')
    AX[0].set_ylabel(u'Déplacement horizontal [mm]',size=14)
    AX[1].set_ylabel(u'Tassement [mm]',size=14)


    fig.tight_layout()
    fig.savefig(figpath+'/'+prob+'_disp_evol')
    plt.close(fig)

    return figpath+'/'+prob+'_disp_evol'
