# -*- coding: cp1252 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from zsoil_tools import zsoil_results as zr

pathname = '../'
pblist = ['M1146_2D']
tvect = range(2,15)


for kf,prob in enumerate(pblist):
    figname = prob
    res = zr(pathname,prob)
    res.read_rcf()
    res.read_his()
    
    tsteps = []
    for kt,step in enumerate(res.steps):
##        if step.time in tvect and step.sf==0 and step.conv_status==-1:
        if step.sf==0 and step.conv_status==-1:
            tsteps.append(kt)
            
    res.out_steps = tsteps
    res.read_dat()
    res.read_s00()
    res.read_s00(res_type='reactions')
    res.read_s03()  # trusses
    res.read_s04()  # beams

    blist = []
    for ke in range(res.nBeams):
        if res.beam.mat[ke] in [11,12]:
            blist.append(ke)
    tlist = []
    for ke in range(res.nTrusses):
        if res.truss.mat[ke] in [15]:
            tlist.append(ke)
            

    fig0 = plt.figure(figsize=(12,10))
    fig1 = plt.figure(figsize=(12,10))
    fig2 = plt.figure(figsize=(12,10))
##    fig.text(0.1,0.97,(u'T=%1.2f'%(tx))+' ('+etapes[tx]+')',size=16,horizontalalignment='left')
    AX = []

    M = [[[1e10,0] for ke in blist],
         [[-1e10,0] for ke in blist]]
    T = [[1e10 for ke in blist],
         [-1e10 for ke in blist]]
    N = [-1e10 for ke in tlist]
    
    for kt in tsteps:
        step = res.steps[kt]
        tx = step.time
        
        for kk,ke in enumerate(blist):
            if step.beam.moment[0][ke]<M[0][kk][0]:
                M[0][kk] = [step.beam.moment[0][ke],
                            step.beam.force[1][ke]]
            elif step.beam.moment[0][ke]>M[1][kk][0]:
                M[1][kk] = [step.beam.moment[0][ke],
                            step.beam.force[1][ke]]
            if step.beam.force[1][ke]<T[0][kk]:
                T[0][kk] = step.beam.force[1][ke]
            elif step.beam.force[1][ke]>T[1][kk]:
                T[1][kk] = step.beam.force[1][ke]
        for kk,ke in enumerate(tlist):
            if step.truss.force[ke]>N[kk]:
                N[kk] = step.truss.force[ke]
        
    AX.append(fig0.add_subplot(111))
    AX.append(fig1.add_subplot(111))
    AX.append(fig2.add_subplot(111))
    s0 = 0.03
    s1 = 0.05
    s2 = 0.01
    lims = [[[1e10,-1e10],[1e10,-1e10]] for k in AX]
    vlim = [[1e10,-1e10] for k in AX]
    for kk,ke in enumerate(blist):
        inel = res.beam.inel[ke]
        xx = np.array([res.coords[0][kn-1] for kn in inel])
        yy = np.array([res.coords[1][kn-1] for kn in inel])
        n = np.array([-yy[1]+yy[0],xx[1]-xx[0]])
        l = np.linalg.norm(n)
        n /= l
        v = [M[0][kk][0]-0.5*l*M[0][kk][1],
             M[0][kk][0]+0.5*l*M[0][kk][1],
             M[1][kk][0]+0.5*l*M[1][kk][1],
             M[1][kk][0]-0.5*l*M[1][kk][1]]
        crds = np.array([[xx[0]+n[0]*s0*v[0],
                          xx[1]+n[0]*s0*v[1],
                          xx[1]+n[0]*s0*v[2],
                          xx[0]+n[0]*s0*v[3]],
                         [yy[0]+n[1]*s0*v[0],
                          yy[1]+n[1]*s0*v[1],
                          yy[1]+n[1]*s0*v[2],
                          yy[0]+n[1]*s0*v[3]]]).T
        lims[0][0][0] = min(lims[0][0][0],np.min(crds,0)[0])
        lims[0][0][1] = max(lims[0][0][1],np.max(crds,0)[0])
        lims[0][1][0] = min(lims[0][1][0],np.min(crds,0)[1])
        lims[0][1][1] = max(lims[0][1][1],np.max(crds,0)[1])
        vlim[0][0] = min(vlim[0][0],min(v))
        vlim[0][1] = max(vlim[0][01],max(v))
        patch = Polygon(crds,fc='none',ec='k')
        AX[0].add_patch(patch)
        AX[0].plot(xx,yy,'k',linewidth=2)
    
        crds = np.array([[xx[0]+n[0]*s1*T[0][kk],
                          xx[1]+n[0]*s1*T[0][kk],
                          xx[1]+n[0]*s1*T[1][kk],
                          xx[0]+n[0]*s1*T[1][kk]],
                         [yy[0]+n[1]*s1*T[0][kk],
                          yy[1]+n[1]*s1*T[0][kk],
                          yy[1]+n[1]*s1*T[1][kk],
                          yy[0]+n[1]*s1*T[1][kk]]]).T
        lims[1][0][0] = min(lims[1][0][0],np.min(crds,0)[0])
        lims[1][0][1] = max(lims[1][0][1],np.max(crds,0)[0])
        lims[1][1][0] = min(lims[1][1][0],np.min(crds,0)[1])
        lims[1][1][1] = max(lims[1][1][1],np.max(crds,0)[1])
        vlim[1][0] = min(vlim[1][0],T[0][kk])
        vlim[1][1] = max(vlim[1][1],T[1][kk])
        patch = Polygon(crds,fc='none',ec='k')
        AX[1].add_patch(patch)
        AX[1].plot(xx,yy,'k',linewidth=2)

    for kk,ke in enumerate(tlist):
        inel = res.truss.inel[ke]
        xx = np.array([res.coords[0][kn-1] for kn in inel])
        yy = np.array([res.coords[1][kn-1] for kn in inel])
        n = np.array([-yy[1]+yy[0],xx[1]-xx[0]])
        l = np.linalg.norm(n)
        n /= l    
        crds = np.array([[xx[0],
                          xx[1],
                          xx[1]+n[0]*s2*N[kk],
                          xx[0]+n[0]*s2*N[kk]],
                         [yy[0],
                          yy[1],
                          yy[1]+n[1]*s2*N[kk],
                          yy[0]+n[1]*s2*N[kk]]]).T
        lims[2][0][0] = min(lims[2][0][0],np.min(crds,0)[0])
        lims[2][0][1] = max(lims[2][0][1],np.max(crds,0)[0])
        lims[2][1][0] = min(lims[2][1][0],np.min(crds,0)[1])
        lims[2][1][1] = max(lims[2][1][1],np.max(crds,0)[1])
        vlim[2][0] = min(vlim[2][0],N[kk])
        vlim[2][1] = max(vlim[2][01],N[kk])
        patch = Polygon(crds,fc='none',ec='k')
        AX[2].add_patch(patch)
        AX[2].plot(xx,yy,'k',linewidth=2)
        
##        ax.set_xlabel(u'Dépl. horiz. [mm]',size=14)
##        
##        ax = fig.add_subplot(2,3,2)
##        ax.plot(MM_g,yb_g)
##        ax.plot([0,0],[yg0,yg1],'k',linewidth=4)
##        AX.append(ax)
##        ax.set_xlabel(u'Moment [kNm/m]',size=14)
##
##        ax = fig.add_subplot(2,3,3)
##        ax.plot(VV_g,yb_g)
##        ax.plot([0,0],[yg0,yg1],'k',linewidth=4)
####        ax.arrow(0,ygt0,step.truss.force[t0],0,head_starts_at_zero=True,head_width=10)
####        ax.arrow(0,ygt1,step.truss.force[t1],0,head_starts_at_zero=True)
##        if tx<8:
##            ax.annotate('  %1.0f kN/m'%(step.truss.force[t0]),xy=(0,ygt0))
##        else:
##            ax.annotate('  %1.0f kN/m'%(sum([step.nodal.r_disp[0][kk] for kk in kgb0])),xy=(0,ygb0))
##        if tx<7:
##            ax.annotate('  %1.0f kN/m'%(step.truss.force[t1]),xy=(0,ygt1))
##        else:
##            ax.annotate('  %1.0f kN/m'%(sum([step.nodal.r_disp[0][kk] for kk in kgb1])),xy=(0,ygb1))
##            ax.annotate('  %1.0f kN/m'%(sum([step.nodal.r_disp[0][kk] for kk in kgb2])),xy=(0,ygb2))
##        AX.append(ax)
##        ax.set_xlabel(u'Effort tranchant [kN/m]',size=14)
##        
##        ax = fig.add_subplot(2,3,4)
##        ax.plot(DD_d,y_d)
##        ax.plot([0,0],[yd0,yd1],'k',linewidth=4)
##        AX.append(ax)
##        ax.set_xlabel(u'Dépl. horiz. [mm]',size=14)
##        
##        ax = fig.add_subplot(2,3,5)
##        ax.plot(MM_d,yb_d)
##        ax.plot([0,0],[yd0,yd1],'k',linewidth=4)
##        AX.append(ax)
##        ax.set_xlabel(u'Moment [kNm/m]',size=14)
##        
##        ax = fig.add_subplot(2,3,6)
##        ax.plot(VV_d,yb_d)
##        ax.plot([0,0],[yd0,yd1],'k',linewidth=4)
##        if tx<8:
##            ax.annotate('  %1.0f kN/m'%(-step.truss.force[t0]),xy=(0,ydt0))
##        else:
##            ax.annotate('  %1.0f kN/m'%(sum([step.nodal.r_disp[0][kk] for kk in kdb0])),xy=(0,ydb0))
##        if tx<7:
##            ax.annotate('  %1.0f kN/m'%(-step.truss.force[t1]),xy=(0,ydt1))
##        else:
##            ax.annotate('  %1.0f kN/m'%(sum([step.nodal.r_disp[0][kk] for kk in kdb1])),xy=(0,ydb1))
##        AX.append(ax)
##        ax.set_xlabel(u'Effort tranchant [kN/m]',size=14)
##        
    for kax,ax in enumerate(AX):
##        ax.set_ylabel(u'Elévation [m]',size=14)
        ax.grid(b=True,which='both',axis='x')
        ax.grid(b=True,which='both',axis='y')
        ax.set_xlim(lims[kax][0])
        ax.set_ylim(lims[kax][01])

    AX[0].annotate('$M_{min}=%1.0f$ kNm/m\n$M_{max}=%1.0f$ kNm/m'%(vlim[0][0],vlim[0][1]),
                   xy=(0.5,0.5),xycoords='axes fraction')
    AX[1].annotate('$T_{min}=%1.0f$ kN/m\n$T_{max}=%1.0f$ kN/m'%(vlim[01][0],vlim[1][1]),
                   xy=(0.5,0.5),xycoords='axes fraction')
    AX[2].annotate('$N_{min}=%1.0f$ kN/m\n$N_{max}=%1.0f$ kN/m'%(vlim[2][0],vlim[2][1]),
                   xy=(0.5,0.5),xycoords='axes fraction')
    AX[2].set_aspect('equal')

    fig0.tight_layout()
    fig0.savefig('paroi_M')
    fig1.tight_layout()
    fig1.savefig('paroi_T')
    fig2.tight_layout()
    fig2.savefig('paroi_N')







