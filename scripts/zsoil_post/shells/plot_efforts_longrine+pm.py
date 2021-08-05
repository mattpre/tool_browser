# -*- coding: cp1252 -*-
# @description Plotting struts with normal forces or section data.
# @input zsoil results
# @output png 
# @author Matthias Preisig
# @date 2018/31/05
import os,math,cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import cPickle as pickle


from zsoil_tools import zsoil_results as zr


pathname = '..'
pblist = ['M1119_3D_v4_ROR_P02']
##pblist = ['M1119_3D_v2_1_offset_mem3']


for kf,prob in enumerate(pblist):
    figname = prob
    try:
        res = pickle.load(open(prob+'_long+pm.p', "rb" ))
        tsteps = res.out_steps
        tvect = [res.steps[kt].time for kt in tsteps]
    except:
        res = zr(pathname,prob)
        res.read_rcf()
        res.read_his()

        tvect = [8]
        tsteps = []
        tsteps_plot = []
        for kt,step in enumerate(res.steps):
            if step.time in tvect and step.sf==0 and step.conv_status==-1:
                tsteps.append(kt)
                
        res.out_steps = tsteps
        res.read_dat()
        res.read_s02()
        res.read_s04()
        pickle.dump(res, open(prob+'_long+pm.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    ypm = [[540,536.5],
           [536.5,534.016667],
           [534.016667,531.216667],
           [531.216667,528.35],
           [528.35,525.35]]
    
    spos = set()
    for ke in range(res.nBeams):
        if res.shell.mat[ke]==6:
            inel = res.shell.inel[ke]
            xx = [res.coords[0][kn-1] for kn in inel]
            yy = [res.coords[1][kn-1] for kn in inel]
            zz = [res.coords[2][kn-1] for kn in inel]
            sp = int('%1.0f'%(np.mean(xx)*10+1000*np.mean(zz)))
            spos.add(sp)
    spos = list(spos)
    
    slists = [[[] for kk in spos] for k in ypm]
    for ke in range(res.nShells):
        if res.shell.mat[ke]==6:
            inel = res.shell.inel[ke]
            xx = [res.coords[0][kn-1] for kn in inel]
            yy = [res.coords[1][kn-1] for kn in inel]
            zz = [res.coords[2][kn-1] for kn in inel]
            sp = int('%1.0f'%(np.mean(xx)*10+1000*np.mean(zz)))
            if sp in spos:
                spind = spos.index(sp)
                ym = np.mean(yy)
                kspr = -1
                for ky,yp in enumerate(ypm):
                    if ym>yp[1] and ym<yp[0]:
                        kspr = ky
                        break
                if kspr>-1:
                    slists[kspr][spind].append(ke)

    ylev = set()
    for ke in range(res.nBeams):
##        if res.beam.mat[ke]==12:
        inel = res.beam.inel[ke]
        yy = [res.coords[1][kn-1] for kn in inel]
        ylev.add(yy[0])
    ylev = list(ylev)
    ylev.sort(reverse=True)
        
    
    ebem = [[0 for kk in spos] for k in ylev]
    nbem = [set() for k in ylev]
    spos2 = set()
    spos3 = [0 for k in spos]
    for ke in range(res.nBeams):
##        if res.beam.mat[ke]==12:
        inel = res.beam.inel[ke]
        xx = [res.coords[0][kn-1] for kn in inel]
        yy = [res.coords[1][kn-1] for kn in inel]
        zz = [res.coords[2][kn-1] for kn in inel]
        sp = int('%1.0f'%(np.mean(xx)*10+1000*np.mean(zz)))
        spos2.add(sp)
        spind = spos.index(sp)
        spos3[spind] += 1
        kspr = ylev.index(np.mean(yy))
        ebem[kspr][spind] = ke
    spos2 = list(spos)

    N = [[0 for ke in range(len(spos))] for kl in range(len(ylev))]
    M = [[0 for ke in range(len(spos))] for kl in range(len(ylev))]
    T = [[0 for ke in range(len(spos))] for kl in range(len(ylev))]
    Ns = [[0 for ke in range(len(spos))] for kl in range(len(ylev))]
    Ms = [[0 for ke in range(len(spos))] for kl in range(len(ylev))]
    Ts = [[0 for ke in range(len(spos))] for kl in range(len(ylev))]

    kt = tsteps[-1]
    step = res.steps[kt]
    
    for kl in range(len(ylev)):
##        for kke,ke in enumerate(etrs[kl]):
##            N[kl][kke] = step.truss.force[ke]
        for kke,ke in enumerate(ebem[kl]):
            N[kl][kke] = step.beam.force[0][ke]
            T[kl][kke] = step.beam.force[1][ke]
            M[kl][kke] = step.beam.moment[2][ke]
        for kke,splist in enumerate(slists[kl]):
            for ke in splist:
                inel = res.shell.inel[ke]
                yy = [res.coords[1][kn-1] for kn in inel]
                dy = max(yy)-min(yy)
                Ns[kl][kke] += step.shell.smforce[1][ke]*dy
                Ts[kl][kke] += step.shell.sqforce[1][ke]*dy
                Ms[kl][kke] += step.shell.smoment[1][ke]*dy

    for kplot in range(3):
        if kplot==0:
            scale = 0.001
            Val = [[N[kl][kp]+Ns[kl][kp] for kp in range(len(spos))] for kl in range(len(ylev))]
##            Val = [[N[kl][kp] for kp in range(len(spos))] for kl in range(len(ylev))]
            lab = 'Normalkraft'
            lab2 = 'N'
            unit = 'kN'
        elif kplot==1:
            scale = 0.002
            Val = [[M[kl][kp]+Ms[kl][kp] for kp in range(len(spos))] for kl in range(len(ylev))]
            lab = 'Biegemoment'
            lab2 = 'M'
            unit = 'kNm'
        else:
            scale = 0.002
            Val = [[T[kl][kp]+Ts[kl][kp] for kp in range(len(spos))] for kl in range(len(ylev))]
            lab = 'Querkraft'
            lab2 = 'T'
            unit = 'kN'
        fig = plt.figure(figsize=(20,12))
        bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
        fig.text(0.1,0.8,u'%s in\nLongarine und Schlitzwand\nbei Aushubende %s'%(lab,unit),
                 size=18,bbox=bbox_props)
        fig.text(0.01,0.01,prob,size=8)

        stit = [u'Spriesslage %i (%1.2f MüM)'%(k+1,ylev[k]) for k in range(len(ylev))]

        for kl in range(len(ylev)):
            ax = fig.add_subplot(2,3,kl+2)
            ax.axis('off')
            for kke,ke in enumerate(ebem[kl]):
                inel = res.beam.inel[ke]
                p0 = np.array([res.coords[0][inel[0]-1],
                               -res.coords[2][inel[0]-1]])
                p1 = np.array([res.coords[0][inel[1]-1],
                               -res.coords[2][inel[1]-1]])
                n = p1-p0
                l = np.linalg.norm(n)
                n /= l
                n = np.array([n[1],-n[0]])
                val = Val[kl][kke]
                p01 = p0 + n*scale*val
                p11 = p1 + n*scale*val
                ang = 180./math.pi*np.angle(n[0]+n[1]*1j)
                ax.plot([p0[0],p01[0],p11[0],p1[0]],
                        [p0[1],p01[1],p11[1],p1[1]],'b')
                ax.plot([p0[0],p1[0]],
                        [p0[1],p1[1]],'k')
                ax.annotate('%1.0f'%(val),xy=0.5*(p0+p1)+1.2*n,
                            rotation=ang,ha='center',va='center')
    ##                        bbox=dict(boxstyle='round,pad=0.2',fc='w',ec='k'))

    ##        for kke in range(len(etrs[kl])):
    ##            inel = res.truss.inel[etrs[kl][kke]]
    ##            X = [res.coords[0][kn-1] for kn in inel]
    ##            Z = [-res.coords[2][kn-1] for kn in inel]
    ##            if X[1]-X[0]==0:
    ##                ang = 90
    ##            else:
    ##                ang = math.atan((Z[1]-Z[0])/(X[1]-X[0]))#/math.pi*180
    ##            ax.plot(X,Z,'k')
    ##            ax.annotate('%1.0f'%(N[kl][kke]),xy=(np.mean(X),np.mean(Z)),
    ##                        rotation=0,ha='center',va='center',
    ##                        bbox=dict(boxstyle='round,pad=0.2',fc='w',ec='k'))
                    
            ax.set_title(stit[kl],
                         size=16)
            xl = ax.get_xlim()
            yl = ax.get_ylim()
            ax.set_xlim(xl[0]-1,xl[1]+1)
            ax.set_ylim(yl[0]-1,yl[1]+1)
    ##        ax.set_xlabel('z [m]',size=16)
    ##        ax.set_ylabel('x [m]',size=16)
            ax.set_aspect('equal')

        fig.tight_layout()
        fig.savefig(prob+'_%slong+pm'%(lab2))
        plt.close(fig)


        

            




