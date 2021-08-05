# -*- coding: cp1252 -*-
# @description Plotting evolution of nodal displacements.
# @output displacement plot
# @author Matthias Preisig
# @date 2017/08/18
import os,math
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
from matplotlib.patches import Polygon
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection


from zsoil_tools import zsoil_results as zr

pathname = '../profil001'
pblist = ['M1161_P001_v2']

etapes = [(2,u'tunnel LO excavé'),
          (3,u'Palace construit'),
          (5,u'tunnel Ouest excavé'),
          (7,u'tunnel Est excavé')]

for kp,prob in enumerate(pblist):

    figname = prob
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
        tsteps_plot = []
        for kt,step in enumerate(res.steps):
            if step.sf==0 and step.conv_status==-1:
                tsteps.append(kt)

        res.out_steps = tsteps
        res.read_dat()
        res.read_s00()
        pickle.dump(res, open(prob+'_disp.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    crds = [np.array([1.74955e+01,2.49552e+01]),
             np.array([1.40613e+01,2.40872e+01]),
             np.array([2.27609e+01,1.73308e+01]),
             np.array([1.22234e+01,1.79991e+01]),
             np.array([1.86667e+01,4.65000e+01]),
             np.array([2.28444e+01,3.15000e+01]),
             np.array([4.00519e+01,3.15000e+01]),
             np.array([-1.81419e+00,2.31425e+01]),
             np.array([2.38765e+01,2.53477e+01])]
    labs = [u'Clé de voûte tunnel LO',
            u'Point de voûte tunnel LO',
            u'Piédroit Est tunnel LO',
            u'Piédroit Ouest tunnel LO',
            u'TN à côté Hôtel Palace, 465 msm',
            u'Base Ouest de fondation Palace',
            u'Base Est de fondation Palace',
            u'Clé de voûte tunnel Ouest',
            u'Clé de voûte tunnel Est']
    bnds = set()
    for ke in range(res.nBeams):
        inel = res.beam.inel[ke]
        for kn in inel:
            bnds.add(kn-1)
    bnds = list(bnds)
    LOnds = set()
    for ke in range(res.nVolumics):
        if res.vol.EF[ke] in [1,12]:
            inel = res.vol.inel[ke]
            for kn in inel:
                LOnds.add(kn-1)
    LOnds = sorted(list(LOnds))
    astr = str()
    for kn in LOnds:
        astr += ' %i'%(kn+1)
    nlist = [-1 for k in labs]
    for kn in range(res.nNodes):
        if not kn in bnds:
            for kk,crd in enumerate(crds):
                if abs(res.coords[0][kn]-crd[0])<1e-3:
                    if abs(res.coords[1][kn]-crd[1])<1e-3:
                        if kk<4:
                            if kn in LOnds:
                                nlist[kk] = kn
                        else:
                            nlist[kk] = kn
    if min(nlist)==-1:
        print('point not found!')
##    step0 = res.steps[tsteps[2]]

    tx = []
    uh = [[] for kn in nlist]
    uv = [[] for kn in nlist]
    for kt in tsteps:
        step = res.steps[kt]
        tx.append(step.time)
        for kk,kn in enumerate(nlist):
            uh[kk].append(1e3*step.nodal.disp[0][kn])
            uv[kk].append(1e3*step.nodal.disp[1][kn])

    lstr = ['bs-','r*-','g<-','mo-']
    lstr = ['b','g','r','m','y','c','b--','g--','r--','m--','y--','c--']

    fig = plt.figure(figsize=(12,12))
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    fig.text(0.01,0.01,prob,size=8)

    axh = fig.add_subplot(211)
    axv = fig.add_subplot(212)
    AX = [axh,axv]

    for kk,u in enumerate(uh):
        axh.plot(tx,u,lstr[kk],label='%i '%(kk+1)+labs[kk])
    for kk,u in enumerate(uv):
        axv.plot(tx,u,lstr[kk],label='%i '%(kk+1)+labs[kk])
        
    axh.set_ylabel(u'Déplacement horizontal [mm]',size=14)
    axv.set_ylabel(u'Déplacement vertical [mm]',size=14)
    axh.legend(loc='best',fontsize=10)
##    ax.set_xlim([0,75])
    for ax in AX:
        ax.grid('on')
##        ax.set_xlabel(u'Etape',size=14)
        ax.set_xticks([e[0] for e in etapes])
        ax.set_xticklabels(['T=%i %s'%(e[0],e[1]) for e in etapes])
        plt.setp(ax.xaxis.get_majorticklabels(),rotation=10,ha='right')

    fig.tight_layout()
    fig.savefig(prob+'_disp')
    plt.close(fig)
    
######################################################

    triangles0 = []
    patches = []
    mat = []
    bounds = [[1e10,-1e10],
              [1e10,-1e10]]
    for ke in range(res.nVolumics):
##        if res.EF[res.vol.EF[ke]][0]<step.time and res.EF[res.vol.EF[ke]][1]>=step.time:
        inel = res.vol.inel[ke]
        triangles0.append([inel[kk]-1 for kk in [0,1,2]])
        triangles0.append([inel[kk]-1 for kk in [0,2,3]])
        xy = np.zeros(4*2)
        xy.shape = (4,2)
        for kk,kn in enumerate(inel):
            xy[kk][0] = res.coords[0][kn-1]
            xy[kk][1] = res.coords[1][kn-1]
        bounds[0] = [min(bounds[0][0],min(xy[:,0])),max(bounds[0][1],max(xy[:,0]))]
        bounds[1] = [min(bounds[1][0],min(xy[:,1])),max(bounds[1][1],max(xy[:,1]))]
        mat.append(res.vol.EF[ke])
        poly = Polygon(np.array(xy))
        patches.append(poly)


    fig = plt.figure(figsize=(8,6))

    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
    fig.text(0.01,0.01,prob,size=8)
    AX = []
    ax = fig.add_subplot(111)
    
    cmap = plt.cm.get_cmap('Set3')
    norm = plt.Normalize(min(mat), max(mat))
    pc = PatchCollection(patches,edgecolors=('none',),norm=norm,cmap=cmap,antialiased=False)
##    pc.set_clim(vmin=minima,vmax=maxima)
    pc.set_array(np.array(mat))
    ax.add_collection(pc)

##    clist = list(set(mat))
##    handles = []
##    for col in clist:
##        handles.append(Polygon([(0,0),(10,0),(0,-10)],color=pc.cmap(pc.norm(col)),
##                               label='%s m/s'%(mat_dict[col])))
##
##    plt.legend(handles=handles)

    for kp in range(len(crds)):
        ax.plot(crds[kp][0],crds[kp][1],'k.')
        ax.annotate('%i'%(kp+1),xy=crds[kp],xytext=(2,2),
                    textcoords='offset points',size=14)

    bounds[1][1] += 5
##    ax.set_xlim(bounds[0])
    ax.set_xlim([-10,45])
    ax.set_ylim(bounds[1])
    ax.set_aspect('equal')
    ax.axis('off')

    bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
    fig.text(0.01,0.01,prob,size=8)

    fig.tight_layout()
    fig.savefig(prob+'_mat')
    plt.close(fig)

    of = open(prob+'_disp.csv','w')
    ktOuest0 = 4
    ktOuest1 = 7
    ktEst0 = 8
    ktEst1 = 11
    tOuest0 = tx[ktOuest0]
    tOuest1 = tx[ktOuest1]
    tEst0 = tx[ktEst0]
    tEst1 = tx[ktEst1]
    of.write('Pt nr;label;x;y;Duh(%1.2g-%1.2f);Duv(%1.2g-%1.2f);Duh(%1.2g-%1.2f);Duv(%1.2g-%1.2f)\n'%(tOuest1,tOuest0,tOuest1,tOuest0,tEst1,tEst0,tEst1,tEst0))
    for kp in range(len(crds)):
        of.write('%i;%s;%1.2f;%1.2f;%1.8e;%1.8e;%1.8e;%1.8e\n'%(kp+1,'',crds[kp][0],crds[kp][1],
                                                    uh[kp][ktOuest1]-uh[kp][ktOuest0],
                                                    uv[kp][ktOuest1]-uv[kp][ktOuest0],
                                                    uh[kp][ktEst1]-uh[kp][ktEst0],
                                                    uv[kp][ktEst1]-uv[kp][ktEst0]))
    of.close()


        

            




