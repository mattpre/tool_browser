# -*- coding: cp1252 -*-
import numpy as np
import math
import matplotlib.pyplot as plt
import cPickle as pickle


from zsoil_tools import zsoil_results as zr

pathname = '..'
pblist = ['M1133_full_v1_4','M1133_full_v1_4_l=100%','M1133_full_v2_l=100%']
pblist = ['M1133_Full_v2_l=100%_E10_VoussMet']
for prob in pblist:
    try:
        res = pickle.load(open(prob+'_disp.p', "rb" ))
        print prob+' loaded'
        tsteps = res.out_steps
        tvect = [res.steps[kt].time for kt in tsteps]
    except:

        res = zr(pathname,prob)
        res.read_rcf()
        res.read_his()
        tsteps = []
        for kt,step in enumerate(res.steps):
            if step.conv_status==-1:
                if abs(step.time%1)<1e-3:
                    tsteps.append(kt)
        res.out_steps = tsteps
        res.read_dat()
##        res.read_LTF()
        res.read_s00()

    shnds = set()
    for ke in range(res.nShells):
        if res.shell.mat[ke]==11:
            inel = res.shell.inel[ke]
            for kn in inel:
                shnds.add(kn-1)
    shnds = list(shnds)

    for ylevel in [22.3,34.87,44.0,12.84,-3.0]:

        nlist = []
        for kn in range(res.nNodes):
            crd = np.array([res.coords[0][kn],res.coords[1][kn],res.coords[2][kn]])
            if abs(crd[1]-ylevel)<1e-3:
                if np.linalg.norm([crd[0],crd[2]])<33.1:
                    if kn in shnds:
                        nlist.append((kn,crd,np.angle(crd[0]-crd[2]*1j)))
        nlist.sort(key=lambda v:v[2])
        nlist = [(v[0],v[1]) for v in nlist]
        nlist.append(nlist[0])
        xn = [v[1][0] for v in nlist]
        zn = [-v[1][2] for v in nlist]
        
        fig = plt.figure(figsize=(12,8))
        fig.text(0.01,0.01,prob,size=8)
        ax = plt.subplot(111)

        etapes = [#u'réalisation pm',
    ##              u'percement Est L16',
    ##              u'percement Ouest L16',
                  u'excavation S1',
    ##              u'percement Est L17N',
    ##              u'percement Ouest L17N',
                  u'excavation S2',
                  u'excavation f.f.']
    ##              u'percement Est L17S',
    ##              u'percement Ouest L17S']
        ktvect = [2,12,14,25,35,37,48,51,62,64]
        ktvect = [25,48,51]

        for step in res.steps:
            if abs(step.time-2)<1e-3:
                step0 = step
                break
        
        scale = 1.e3
        for kkt,kt in enumerate(res.out_steps):
            step = res.steps[kt]
            if step.time in ktvect:
                ind = ktvect.index(step.time)
                dxn = [step.nodal.disp[0][v[0]]-step0.nodal.disp[0][v[0]] for v in nlist]
                dzn = [-step.nodal.disp[2][v[0]]+step0.nodal.disp[2][v[0]] for v in nlist]
                ax.plot([xn[kk]+scale*dxn[kk] for kk in range(len(nlist))],
                        [zn[kk]+scale*dzn[kk] for kk in range(len(nlist))],
                        label=etapes[ind])
                
    ##        ax.grid('on')
        ax.set_aspect('equal')
        ax.legend()
        
        for ke in range(res.nShells):
            if res.shell.mat[ke]==11:
                inel = res.shell.inel[ke]
                xx = [res.coords[0][kn-1] for kn in inel]
                yy = [res.coords[1][kn-1] for kn in inel]
                zz = [-res.coords[2][kn-1] for kn in inel]
                if max(yy)>43.99:
                    ax.plot(xx,zz,'k')

        R = 33.07
        for kk in range(1,6):
            xx = [(R-kk*scale*2e-3)*math.cos(ktheta/999.*2*math.pi) for ktheta in range(1000)]
            zz = [(R-kk*scale*2e-3)*math.sin(ktheta/999.*2*math.pi) for ktheta in range(1000)]
            ax.plot(xx,zz,'k--',linewidth=0.5)
            ax.annotate('%i mm'%(kk*2),xy=(xx[0],zz[0]),rotation=90,ha='right')
    ##        ax.set_xlim([0,4])
    ##        ax.set_ylim([-1.5,5.5])
            
        ax.annotate(u'Déformée de la paroi moulée à %1.1f NGF'%(ylevel),
                    xy=(0,0),ha='center',size=12)
        
        fig.tight_layout()
        fig.savefig(prob+'_disp_pm_%1.0fNGF'%(ylevel))
        plt.close(fig)
