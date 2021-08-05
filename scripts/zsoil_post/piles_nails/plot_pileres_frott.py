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


pathname = '../..'
pblist = [#'M1013_O2021_v1_pilesV1017_cent1_L1Ost']
##          'M1013_O2021_v1_pilesV1017_cent1_L2Ost',
          'M1013_O2021_v1_piles_red1_L1Ost']
##          'M1013_O_v8_pilesV1017_cent2_l171110']


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

        tvect = [4,5]
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
##    pilecrds = [(c[1],-c[0]) for c in pilecrds]
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

    FNsc = 5e-4
    diam = {13:0.9,14:1.0}
    fig = plt.figure(figsize=(8,12))
    fig.text(0.01,0.01,prob,size=8)
    v0 = np.array([0,0])
    v1 = np.array([-1,0.5])
    v2 = np.array([-1,1.5])
    ax = fig.add_subplot(111)
    for k in range(len(pilecrds)):
        pos = pilecrds[k]
        ax.add_patch(Rectangle(pos+v0*0,
                               width=0.5,height=-FN[k][0]*FNsc,color=colvect[0]))
        ax.annotate('%1.0f kN'%(FN[k][0]),xy=pos+v0,
                    ha='right',rotation=90,size=7)

        mat = res.beam.mat[elist[k][0]]
        qs = 70
        qp = -2000
        cel = clist[k]

        L = 0
        qratio = 0
##        print(pilecrds[k])
        for kke,ke in enumerate(cel):#[1:]):
            vke = res.piles[k].volumics[kke]
            vmat = res.vol.mat[vke[0]]
            if vmat==1:
                qs = 0
            elif vmat==2:
                qs = 170
            elif vmat==3:
                qs = 110
            elif vmat==5:
                qs = 65
            else:
                print('pile in layer %i'%(vmat))
            
            dl = abs(res.coords[1][res.cnt.inel[ke][0]-1]-res.coords[1][res.cnt.inel[ke][1]-1])
            L += dl
            qratio += max(0,np.mean([step.cnt.stress[0][ke][kgp]/qs for kgp in [0,1]]))*dl
##            if abs(pilecrds[k][0]-66)<0.1:
##                print(min(1,max(0,np.mean([step.cnt.stress[0][ke][kgp]/qs for kgp in [0,1]]))),mat)
##            qratio += np.mean([step.cnt.str_level[ke][kgp] for kgp in [0,1]])*dl
        vke = res.piles[k].volumics[-1]
        vmat = res.vol.mat[vke[0]]
        if vmat==2:
            qp = -4000
        elif vmat==3:
            qp = -2500
        elif vmat==5:
            qp = -1000
        else:
            print('tip in layer %i'%(vmat))
        p = ax.pie([qratio/L,1-qratio/L],center=pos+v1,
                   radius=0.5,colors=[colvect[1],'grey'])
        p[0][1].set_alpha(0.4)
        p = ax.pie([step.cnt.stress[2][ctlist[k]][0]/qp,1-step.cnt.stress[2][ctlist[k]][0]/qp],
                   center=pos+v2,radius=0.5,colors=[colvect[2],'grey'])
        p[0][1].set_alpha(0.4)

    p0 = np.array([200,-355])
    p1 = np.array([202,-355])
    ax.add_patch(Rectangle(p0+v0,width=0.5,height=500*FNsc,color=colvect[0]))
    ax.annotate('%1.0f kN'%(-500),xy=p0-v0*0.2,
                ha='right',rotation=90,size=7)
    p = ax.pie([0.25,0.75],center=p0+v1,radius=0.5,colors=[colvect[1],'grey'])
    p[0][1].set_alpha(0.4)
    p = ax.pie([0.20,0.80],center=p0+v2,radius=0.5,colors=[colvect[2],'grey'])
    p[0][1].set_alpha(0.4)
    ax.annotate('Normalkraft Pfahlkopf',xy=p1+v0,va='center',size=8)
    ax.annotate('Mobilisierte Mantelreibung (Annahme: $q_s$=70 kPa)',xy=p1+v1,va='center',size=8)
    ax.annotate('Mobilisierter Spitzenwiderstand (Annahme: $q_p$=2500 kPa)',xy=p1+v2,va='center',size=8)

    ax.set_xlim(bounds[0])
    ax.set_ylim(bounds[1])
    fig.tight_layout()
    fig.savefig(prob+'_pileres_frott')
    plt.close(fig)








