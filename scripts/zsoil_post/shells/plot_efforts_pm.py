# -*- coding: cp1252 -*-
import numpy as np
import math
import cPickle as pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from zsoil_tools import zsoil_results as zr

pathname = '..'
pblist = ['M1133_Full_v1_4','M1133_full_v1_4_l=100%','M1133_full_v2_l=100%']
pblist = ['M1133_Full_v2_l=100%_E10_VoussMet']

labstr = ['Effort normal horizontal [kN/m]',
          'Effort normal vertical [kN/m]',
          'Moment transversal [kNm/m]',
          'Moment usuel [kNm/m]',
          'Eff. tranchant transversal [kN/m]',
          'Eff. tranchant usuel [kN/m]']

for kf,prob in enumerate(pblist):
    figname = prob
    try:
        res = pickle.load(open(prob+'_shell.p', "rb" ))
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
        res.read_s02()
        pickle.dump(res, open(prob+'_shell.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    tsteps_plot = []
    for kt in tsteps:
        if res.steps[kt].time in [12,14,25,35,37,48,51,62,64,80]:
            tsteps_plot.append(kt)
    etapes = [#u'réalisation pm',
              u'percement Est L16',
              u'percement Ouest L16',
              u'excavation S1',
              u'percement Est L17N',
              u'percement Ouest L17N',
              u'excavation S2',
              u'excavation f.f.',
              u'percement Est L17S',
              u'percement Ouest L17S',
              u'fin de simulation']

    ymin = 0
    for inel in res.shell.inel:
        for kn in inel:
            ymin = min(ymin,res.coords[1][kn-1])


##if 1==1:
    R = 180/math.pi
    minima = [[1e10 for kk in range(res.nShells)] for k in range(6)]
    maxima = [[-1e10 for kk in range(res.nShells)] for k in range(6)]
    
    for kkt,kt in enumerate(tsteps_plot):
        step = res.steps[kt]
        tx = step.time
        
        fig = plt.figure(figsize=(18,16))
        bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
        fig.text(0.01,0.01,prob,size=8)
        AX = []

        patches = []
        cvect = [[] for k in range(6)]
        for ke in range(res.nShells):
            if res.shell.mat[ke]==11 and res.EF[res.shell.EF[ke]][1]>=step.time:
                xy = np.zeros(4*2)
                xy.shape = (4,2)
                xx = [res.coords[0][kn-1] for kn in res.shell.inel[ke]]
                yy = [res.coords[1][kn-1] for kn in res.shell.inel[ke]]
                zz = [res.coords[2][kn-1] for kn in res.shell.inel[ke]]
                if (np.mean(xx)**2+np.mean(zz)**2)**0.5<33.1:
                    if np.mean(zz)>0:
                        angle = [abs(180./math.pi*np.angle(xx[kk]+zz[kk]*1j)) for kk in range(4)]
                    else:
                        angle = [-abs(180./math.pi*np.angle(xx[kk]+zz[kk]*1j)) for kk in range(4)]
                    xy = np.array([angle,yy]).T
                    cvect[0].append(step.shell.smforce[1][ke])
                    cvect[1].append(step.shell.smforce[0][ke])
                    cvect[2].append(step.shell.smoment[1][ke])
                    cvect[3].append(step.shell.smoment[0][ke])
                    cvect[4].append(step.shell.sqforce[1][ke])
                    cvect[5].append(step.shell.sqforce[0][ke])
                        
                    poly = Polygon(np.array(xy))
                    patches.append(poly)

        for k in range(6):
            minima[k] = [min(minima[k][kk],cvect[k][kk]) for kk in range(len(cvect[0]))]
            maxima[k] = [max(maxima[k][kk],cvect[k][kk]) for kk in range(len(cvect[0]))]
            pc = PatchCollection(patches,edgecolors=('none',))
            pc.set_array(np.array(cvect[k]))
            ax = fig.add_subplot(3,2,k+1)
            ax.add_collection(pc)
            cb = fig.colorbar(pc)
            cb.set_label(labstr[k],size=16)
            cb.ax.tick_params(labelsize=16)
##            ax.plot([-180,180],[y0,y0],'k')
##            ax.plot([-dload,-dload,dload,dload,-dload],[y1,y0,y0,y1,y1],'k')
            ax.text(-170,-2,'T=%1.0f, %s'%(tx,etapes[kkt]),size=14)
            AX.append(ax)
        
        for ax in AX:
            ax.plot([-180,180],[36.15,36.15],'k')
            ax.plot([-180,180],[28.47,28.47],'k')
            ax.plot([-180,180],[16.8,16.8],'k')
            ax.annotate('L16 Est',xy=(0,22.3),rotation=90,ha='center',va='center')
            ax.annotate('L17S Est',xy=(37.42,22.3),rotation=90,ha='center',va='center')
            ax.annotate('L17N Est',xy=(-37.42,22.3),rotation=90,ha='center',va='center')
            ax.annotate('L17S Ouest',xy=(142.58,22.3),rotation=90,ha='center',va='center')
            ax.annotate('L17N Ouest',xy=(-142.58,22.3),rotation=90,ha='center',va='center')
            ax.set_xlim(-180,180)
            ax.set_ylim(ymin,44)
            ax.set_xlabel(u'Angle p.r. à l\'axe L16 [°]',size=16)
            ax.set_ylabel('Altitude [NGF]',size=16)
    ##        ax.grid(b=True,which='both',axis='x')
    ##        ax.grid(b=True,which='major',axis='y')

        fig.tight_layout()
        fig.savefig(figname+'_carte_t%i_%i'%(tx,(tx%1)*10))
        plt.close(fig)

##    # enveloppes:
##    fig = plt.figure(figsize=(16,16))
##    for k in range(6):
##        pc = PatchCollection(patches,edgecolors=('none',))
##        pc.set_array(numpy.array(minima[k]))
##        ax = fig.add_subplot(3,2,k+1)
##        ax.add_collection(pc)
##        cb = fig.colorbar(pc)
##        cb.set_label(labstr[k],size=16)
##        cb.ax.tick_params(labelsize=16)
##        ax.plot([-R*math.pi,R*math.pi],[-14.85,-14.85],'k')
##        ax.plot([-20.77,-20.77,20.77,20.77,-20.77],[-9.75,-15.35,-15.35,-9.75,-9.75],'k')
##        ax.text(-170,-2,'min')
##        AX.append(ax)
##    
##    for ax in AX:
##        ax.set_xlim(-R*math.pi,R*math.pi)
##        ax.set_ylim(ymin,0)
##        ax.set_xlabel(u'Winkel von Angriffspunkt Presskraft [°]',size=16)
##        ax.set_ylabel('Tiefe [m]',size=16)
##
##    fig.tight_layout()
##    fig.savefig('shell_carte_'+figname+'_min')
##
##    fig = plt.figure(figsize=(16,16))
##    for k in range(6):
##        pc = PatchCollection(patches,edgecolors=('none',))
##        pc.set_array(numpy.array(maxima[k]))
##        ax = fig.add_subplot(3,2,k+1)
##        ax.add_collection(pc)
##        cb = fig.colorbar(pc)
##        cb.set_label(labstr[k],size=16)
##        cb.ax.tick_params(labelsize=16)
##        ax.plot([-R*math.pi,R*math.pi],[-14.85,-14.85],'k')
##        ax.plot([-20.77,-20.77,20.77,20.77,-20.77],[-9.75,-15.35,-15.35,-9.75,-9.75],'k')
##        ax.text(-170,-2,'max')
##        AX.append(ax)
##    
##    for ax in AX:
##        ax.set_xlim(-R*math.pi,R*math.pi)
##        ax.set_ylim(ymin,0)
##        ax.set_xlabel(u'Winkel von Angriffspunkt Presskraft [°]',size=16)
##        ax.set_ylabel('Tiefe [m]',size=16)
##
##    fig.tight_layout()
##    fig.savefig('shell_carte_'+figname+'_max')





