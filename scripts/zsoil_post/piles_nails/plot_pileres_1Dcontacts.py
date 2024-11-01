# -*- coding: cp1252 -*-
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pickle
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
colvect = ['#1f77b4','#ff7f0e','#2ca02c']

from zsoil_tools import zsoil_results as zr


pathname = '..'
pblist = ['M1279_Scheibe_v4b_D1']


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
        pickle.dump(res, open(prob+'.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    cntnodes = [set(),set()]    #[volumic nodes, beam nodes]
    for ke in range(res.nContacts):
        if res.cnt.type[ke]==2:
            cntnodes[0].add(res.cnt.inel[ke][0])
            cntnodes[0].add(res.cnt.inel[ke][1])
            cntnodes[1].add(res.cnt.inel[ke][2])
            cntnodes[1].add(res.cnt.inel[ke][3])
    cntnodes = [list(cntnodes[0]),list(cntnodes[1])]
    pilecrds = []
    plist = []
    for ke in range(res.nBeams):
        inel = res.beam.inel[ke]
        # check if beam has pile interface:
        if inel[0] in cntnodes[1] and inel[1] in cntnodes[1]:
            plist.append(ke)
            # hypothèse: all piles are vertical along y
            pcrd = (res.coords[0][inel[0]-1],res.coords[2][inel[0]-1])
            tag = 0
            for pcrd1 in pilecrds:
                if abs(pcrd[0]-pcrd1[0])+abs(pcrd[1]-pcrd1[1])<1e-3:
                    tag = 1
                    break
            if tag==0:
                pilecrds.append(pcrd)
                    

    elist = []
    clist = []
    ctlist = []
    for kp in range(len(pilecrds)):
        pcrd = pilecrds[kp]
        el = []
        cl = []
        ctl = -1
        for ke in plist:
            inel = res.beam.inel[ke]
            x = res.coords[0][inel[0]-1]
            z = res.coords[2][inel[0]-1]
            if abs(pcrd[0]-x)+abs(pcrd[1]-z)<1e-3:
                el.append(ke)
        for ke in range(res.nContacts):
            inel = res.cnt.inel[ke]
            x = res.coords[0][inel[0]-1]
            z = res.coords[2][inel[0]-1]
            if res.cnt.type[ke]==2:
                if abs(pcrd[0]-x)+abs(pcrd[1]-z)<1e-3:
                    cl.append(ke)
            elif res.cnt.type[ke]==3:
                if abs(pcrd[0]-x)+abs(pcrd[1]-z)<1e-3:
                    ctl = ke
        elist.append(el)
        clist.append(cl)
        if ctl>-1:
            ctlist.append(ctl)


    kntop = []  # 1st node from top of pile
    for kp in range(len(pilecrds)):
        kntop.append(res.beam.inel[elist[kp][0]][0])
    # get pile cap settlement:
    step0 = res.steps[tsteps[0]]
    step = res.steps[tsteps[1]]
    dy_cap = []
    for kp in range(len(pilecrds)):
        dy_cap.append(min(0,step.nodal.disp[1][kntop[kp]-1]-step0.nodal.disp[1][kntop[kp]-1]))
    

    # get normal forces of piles:
    FN = []
    Fs = []
    Ft = []
    for el in elist:
        fn = []
        for ke in el:
            fn.append(step.beam.force[0][ke])
        FN.append(fn)
    for cl in clist:
        fs = []
        for ke in cl:
            fs.append([step.cnt.stress[0][ke][kgp] for kgp in [0,1]])
        Fs.append(fs)
    for ke in ctlist:
        Ft.append(step.cnt.stress[0][ke])

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


   
# mit Kraftverlauf:
# tourner de 90°:
    pilecrds = [(c[1],-c[0]) for c in pilecrds]
# get bounds:
    bounds = [[1e10,-1e10],[1e10,-1e10]]
    for k in range(len(pilecrds)):
        bounds[0][0] = min(bounds[0][0],pilecrds[k][0])
        bounds[0][1] = max(bounds[0][1],pilecrds[k][0])
        bounds[1][0] = min(bounds[1][0],pilecrds[k][1])
        bounds[1][1] = max(bounds[1][1],pilecrds[k][1])
    bounds[0][0] -= 1
    bounds[0][1] += 1
    bounds[1][0] -= 4
    bounds[1][1] += 2

    FNsc = 1e-3
    diam = {11:0.3,12:0.45,13:0.6}
    fig = plt.figure(figsize=(18,3))
    fig.text(0.01,0.01,prob,size=8)
    ax = fig.add_subplot(111)
    FN_stutzen_max = 0
    FN_feld_max = 0
    for k in range(len(pilecrds)):
        pos = pilecrds[k]
        ax.add_patch(Rectangle(pos+np.array([0,0]),
                               width=0.5,height=-FN[k][0]*FNsc,color=colvect[0]))
        
        mat = res.beam.mat[elist[k][0]]
        if mat==15: # pile on symmetry axis --> Fn*2
            ax.annotate('%1.0f kN'%(FN[k][0]*2),xy=pos+np.array([0,-1]),
                        ha='right',rotation=90,size=7)
            FN_feld_max = min(FN_feld_max,FN[k][0]*2)
        elif pos[1] > 1.2:
            ax.annotate('%1.0f kN'%(FN[k][0]),xy=pos+np.array([0,-1]),
                        ha='right',rotation=90,size=7)
            FN_feld_max = min(FN_feld_max,FN[k][0])
        else:
            ax.annotate('%1.0f kN'%(FN[k][0]),xy=pos+np.array([0,0]),
                        ha='right',rotation=90,size=7)
            FN_stutzen_max = min(FN_stutzen_max,FN[k][0])

        if mat==14 or mat==22:
            qs = 70
            qp = -2500
        else:
            qs = 35
            qp = -1250
        cel = clist[k]

        L = 0
        qratio = 0
##        print(pilecrds[k])
        for kke,ke in enumerate(cel):#[1:]):
            dl = abs(res.coords[1][res.cnt.inel[ke][0]-1]-res.coords[1][res.cnt.inel[ke][1]-1])
            L += dl
            qratio += max(0,np.mean([step.cnt.stress[0][ke][kgp]/qs for kgp in [0,1]]))*dl
##            if abs(pilecrds[k][0]-66)<0.1:
##                print(min(1,max(0,np.mean([step.cnt.stress[0][ke][kgp]/qs for kgp in [0,1]]))),mat)
##            qratio += np.mean([step.cnt.str_level[ke][kgp] for kgp in [0,1]])*dl
        p = ax.pie([qratio/L,1-qratio/L],center=pos+np.array([0,-0.5]),
                   radius=0.5,colors=[colvect[1],'grey'])
        p[0][1].set_alpha(0.4)
        if ctlist:
            p = ax.pie([step.cnt.stress[2][ctlist[k]][0]/qp,1-step.cnt.stress[2][ctlist[k]][0]/qp],
                       center=pos+np.array([0,-1.5]),radius=0.5,colors=[colvect[2],'grey'])
        else:
            p = ax.pie([0,1],
                       center=pos+np.array([0,-1.5]),radius=0.5,colors=[colvect[2],'grey'])
        p[0][1].set_alpha(0.4)

    ax.add_patch(Rectangle((33,-4),width=0.5,height=500*FNsc,color=colvect[0]))
    ax.annotate('%1.0f kN'%(-500),xy=(33,-4),
                ha='right',rotation=90,size=7)
    p = ax.pie([0.25,0.75],center=(60,-3.5),radius=0.5,colors=[colvect[1],'grey'])
    p[0][1].set_alpha(0.4)
    p = ax.pie([0.20,0.80],center=(90,-3.5),radius=0.5,colors=[colvect[2],'grey'])
    p[0][1].set_alpha(0.4)
    ax.annotate('Normalkraft Pfahlkopf (x2 auf Symmetrieachse) \n Min Fn - Stützen %1.0f kN - Feld %1.0f kN'%(FN_stutzen_max,FN_feld_max),xy=(34,-3.5),va='center',size=8)
    # ax.annotate('Min Fn - Stutzen - %1.0f kN - Feld - %1.0f kN'%(FN_stutzen_max,FN_feld_max),xy=(34,-4.5),va='center',size=8)
    ax.annotate('Mobilisierte Mantelreibung (Annahme: $q_s$=70 kPa)',xy=(61,-3.5),va='center',size=8)
    ax.annotate('Mobilisierter Spitzenwiderstand (Annahme: $q_p$=2500 kPa)',xy=(91,-3.5),va='center',size=8)
    print(FN_stutzen_max,FN_feld_max)

    ax.set_xlim(bounds[0])
    ax.set_ylim(bounds[1])
    fig.tight_layout()
    fig.savefig(prob+'_pileres_frott')
    plt.close(fig)
