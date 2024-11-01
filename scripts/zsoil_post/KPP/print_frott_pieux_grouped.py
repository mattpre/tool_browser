# -*- coding: cp1252 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pickle
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

from zsoil_tools import zsoil_results as zr

pathname = '../..'
modele = 'M1013_O_2024_v0_5_v24'
pblist = [#modele,
##          modele+'_mitGelb',
          modele+'_piles+10']
##          modele+'_piles+10_DN900',
##          modele+'_pilesv2+10',
##          modele+'_pilesv1_2']
##          modele+'_qsmin_piles+10']
##          modele+'_qsmax_piles+10']


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

        tvect = [3,6]
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

    step0 = res.steps[tsteps[0]]
    step = res.steps[tsteps[1]]

    lengths = set()
    pilecrds = []
    elist = []
    clist = []
    ctlist = []
    dy_cap = []
    volmat = []
    for kp,pile in enumerate(res.piles):
        blist = [res.num_beams.index(ke) for ke in pile.beams]
        elist.append(blist)
        
        yy = [[res.coords[1][kn-1] for kn in res.beam.inel[ke]] for ke in blist]
        lengths.add(max(max(yy))-min(min(yy)))

        xx = [res.coords[0][kn-1] for kn in res.beam.inel[blist[0]]]
        zz = [res.coords[2][kn-1] for kn in res.beam.inel[blist[0]]]
        pilecrds.append(np.array([np.mean(xx),np.mean(zz)]))

        clist.append([res.num_contacts.index(ke) for ke in pile.cnt1D])
        ctlist.append(res.num_contacts.index(pile.cnt0D))

        kntop = res.beam.inel[elist[kp][0]][0]-1      
        dy_cap.append(min(0,step.nodal.disp[1][kntop]-step0.nodal.disp[1][kntop]))

        volmat.append([res.vol.mat[res.num_volumics.index(ke[0])] for ke in pile.volumics])
    lengths = sorted(list(lengths))

    # get normal forces of piles:
    FN = []
    Fs = []
    Ft = []
    for el in elist:
        fn = []
        for ke in el:
            fn.append(step.beam.force[0][ke])
        FN.append(fn)
    # get sleeve friction:
    for cl in clist:
        fs = []
        for ke in cl:
            fs.append([step.cnt.stress[0][ke][kgp] for kgp in [0,1]])
        Fs.append(fs)
    # get mobilized tip resistance:
    for ke in ctlist:
        Ft.append(step.cnt.stress[0][ke])

   
    # Plot pile results (FN, qs,mob, slip) for piles grouped according to length:
    qsdict = {2:160,3:100}
    nCat = len(lengths)
    fig = plt.figure(figsize=(18,14))
    fig.text(0.01,0.01,prob,size=8)
    AX = [[fig.add_subplot(3,nCat,kl*nCat+kk+1) for kl in range(nCat)] for kk in range(3)]
    qpmax = [0 for kk in range(nCat)]
    yqp = [0 for kk in range(nCat)]
    for kpile in range(res.nPiles):
        yy = [[res.coords[1][kn-1] for kn in res.beam.inel[res.num_beams.index(ke)]] for ke in res.piles[kpile].beams]
        length = max(max(yy))-min(min(yy))
        if length<lengths[0]+0.5:
            kl = 0
        elif length<lengths[1]+0.5:
            kl = 1
        else:
            kl = 2
        ax0 = AX[kl][0]
        ax1 = AX[kl][1]
        ax2 = AX[kl][2]
        el = elist[kpile]
        for kke,ke in enumerate(el):
            ax0.plot([FN[kpile][kke],FN[kpile][kke]],
                     [res.coords[1][res.beam.inel[ke][0]-1],
                      res.coords[1][res.beam.inel[ke][1]-1]],'k',alpha=0.2,
                     lw=2,solid_capstyle="butt")
        cel = clist[kpile]
        for kke,ke in enumerate(cel):
            if max([step.cnt.pla_code[ke][kgp] for kgp in [0,1]])==16:
                linecolor = 'r'
            else:
                linecolor = 'b'
            ax1.plot([step.cnt.stress[0][ke][kgp] for kgp in [0,1]],
                     [res.coords[1][res.cnt.inel[ke][0]-1],
                      res.coords[1][res.cnt.inel[ke][1]-1]],linecolor,alpha=0.2,
                     lw=2,solid_capstyle="butt")
            ax2.plot([step.cnt.strain[0][ke][kgp]*1000 for kgp in [0,1]],
                     [res.coords[1][res.cnt.inel[ke][0]-1],
                      res.coords[1][res.cnt.inel[ke][1]-1]],linecolor,alpha=0.2,
                     lw=2,solid_capstyle="butt")
        ax1.barh([369],[-step.cnt.stress[2][ctlist[kpile]][0]*0.1],color='k',
                 alpha=0.2)
        qpmax[kl] = max(qpmax[kl],-step.cnt.stress[2][ctlist[kpile]][0])
        yqp[kl] = min(min(yy))

    for kl in range(3):
        AX[kl][1].annotate('Max. mobilisierter\nSpitzenwiderstand: %1.0f kPa'%(qpmax[kl]),
                           xy=(qpmax[kl]*0.1,367.5),ha='right',va='top')
        AX[kl][0].invert_xaxis()
        AX[kl][0].grid(which='both')
        AX[kl][1].grid(which='both')
        AX[kl][2].grid(which='both')
        AX[kl][0].set_xlabel('Normalkraft [kN]')
        AX[kl][1].set_xlabel(u'Mobil.Mantelreibung [kPa]')
        AX[kl][2].set_xlabel(u'Axiale Relativverschiebung [mm]')
        for kk in range(3):
            AX[kl][kk].set_ylim([398-max(lengths)-2,398])
        AX[kl][1].set_xlim([0,260])
##        AX[kl][2].set_xlim([0,9])
        # Ranges of qs given by geologist:
        AX[kl][1].add_patch(Rectangle(xy=(160,380),width=20,height=20,fc='g',
                                      alpha=0.2,ec=None))
        AX[kl][1].add_patch(Rectangle(xy=(100,350),width=20,height=30,fc='g',
                                      alpha=0.2,ec=None))
    fig.tight_layout()
    fig.savefig(prob+'_frott.svg')
    plt.close(fig)








