# -*- coding: cp1252 -*-
# @description Plotting shell forces envelopes on cut section using vtk.
# @input zsoil results
# @output png
# @author Matthias Preisig
# @date 2018/06/22
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import vtk
import os

from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools


pathname = '..'
pblist = ['M1133_full_v1_4','M1133_full_v1_4_l=100%','M1133_full_v2_l=100%']
pblist = ['M1133_Full_v2_l=100%_E10_VoussMet']
tx = []

for prob in pblist:
    try:
        f = open(pathname+'/pv/'+prob+'.pvd')
        for line in f:
            if 'shell.vtu' in line:
                t = float(line[line.find('timestep=')+9:].split('"')[1])
                tx.append(t)
        print(pathname+'/pv/'+prob+'.pvd loaded')
    except:
        res = zr(pathname,prob)
        res.read_rcf()
        res.read_his()
        tsteps = []
        for kt,step in enumerate(res.steps):
            if step.conv_status in [-1]:
                if abs(step.time-int(step.time))<0.001:
                    tsteps.append(kt)
                    tx.append(step.time)
        res.out_steps = tsteps
        res.read_dat()
        res.read_s00()
        for lab in res.ele_group_labels:
            if lab=='SHELLS':
                res.read_s02()  # shells
        ptemp = os.getcwd()
        os.chdir('../pv')
        vtktools.write_vtu(res,shells=True,verbose=False)
        os.chdir(ptemp)

    zvect = [(30*k)/180.*math.pi for k in range(6)]
    for kz in range(len(zvect)):
        plane = vtk.vtkPlane()
        plane.SetOrigin((0,0,0))
        plane.SetNormal((math.cos(zvect[kz]),0,math.sin(zvect[kz])))

        segposlist = []
        uniquesegments = []
        Tsegments = []

        for kt in range(len(tx)):
            t = tx[kt]
            tstr = '_'+str(int(t)).rjust(3,'0')+'_'+str(int((t-int(t))*100)).rjust(2,'0')
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName('../pv/'+prob+tstr+'_shell.vtu')
            reader.Update()
            mesh = reader.GetOutput()
            
            segments = vtktools.get_section(mesh,plane,EFlist=[1],matlist=[11])
            Tsegments.append(segments)

            for ke,seg in enumerate(segments):
                x = [seg[0][0],seg[1][0]]
                y = [seg[0][1],seg[1][1]]
                segpos = float('%1.2f'%(np.mean(x)+1000*np.mean(y)))
                if segpos in segposlist:
                    kseg = segposlist.index(segpos)
                else:
                    segposlist.append(segpos)
                    uniquesegments.append((x,y))
                    kseg = len(segposlist)-1

            
            if abs(t-2)<0.01:
                # get outline of pm:
                outline_segments = []
                for ke in range(mesh.GetNumberOfCells()):
                    cell = mesh.GetCell(ke)
                    line = []
                    for kk in range(4):
                        pt = mesh.GetPoint(cell.GetPointId(kk))
                        if pt[1]>43.99:
                            line.append(pt)
                    if len(line)==2:
                        outline_segments.append(([line[0][0],line[1][0]],
                                                 [-line[0][2],-line[1][2]]))

        Vals = [[[],[]],[[],[]],[[],[]]]
        for ks in range(len(segposlist)):
            Vals[0][0].append(1e10)
            Vals[0][1].append(-1e10)
            Vals[1][0].append(1e10)
            Vals[1][1].append(-1e10)
            Vals[2][0].append(1e10)
            Vals[2][1].append(-1e10)

        for kt in range(len(tx)):
            segments = Tsegments[kt]
            for kks in range(len(segments)):
                x = [segments[kks][0][0],segments[kks][1][0]]
                y = [segments[kks][0][1],segments[kks][1][1]]
                segpos = float('%1.2f'%(np.mean(x)+1000*np.mean(y)))
                ks = segposlist.index(segpos)
                Vals[0][0][ks] = min(Vals[0][0][ks],segments[kks][2][0][0])
                Vals[0][1][ks] = max(Vals[0][1][ks],segments[kks][2][0][0])
                Vals[1][0][ks] = min(Vals[1][0][ks],segments[kks][2][1][0])
                Vals[1][1][ks] = max(Vals[1][1][ks],segments[kks][2][1][0])
                Vals[2][0][ks] = min(Vals[2][0][ks],segments[kks][2][2][0])
                Vals[2][1][ks] = max(Vals[2][1][ks],segments[kks][2][2][0])
                
        fig = plt.figure(figsize=(10,12))
        bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
        fig.text(0.01,0.01,prob,size=8)


        AX = [fig.add_subplot(221),
              fig.add_subplot(222),
              fig.add_subplot(223),
              fig.add_subplot(224)]

        scale = [-2.e1/max(abs(min(vals[0])),abs(max(vals[1]))) for vals in Vals]
        M0 = []
        M1 = []
        for ks in range(len(segposlist)):
            x = uniquesegments[ks][0]
            y = uniquesegments[ks][1]
            vec = np.array([y[1]-y[0],-x[1]+x[0]])
            norm = np.linalg.norm(vec)
            l2 = 0#-0.5*norm*np.sign(vec[0])
            if norm>1e-3 and abs(y[0]-y[1])>0.01 and abs(abs(np.mean(x))-33.07)<0.05:
                vec /= norm
                if vec[0]>0:
                    ha0 = 'left'
                    ha1 = 'right'
                else:
                    ha1 = 'left'
                    ha0 = 'right'
                for kax,ax in enumerate(AX[1:]):
                    ax.plot(x,y,'k',linewidth=1)
                    if kax==0:
                        v00 = (Vals[0][0][ks]+Vals[2][0][ks]*l2)*vec*scale[kax]
                        v10 = (Vals[0][0][ks]-Vals[2][0][ks]*l2)*vec*scale[kax]
                        M0.append(Vals[0][0][ks]+Vals[2][0][ks]*l2)
                        M0.append(Vals[0][0][ks]-Vals[2][0][ks]*l2)
                        v01 = (Vals[0][1][ks]+Vals[2][1][ks]*l2)*vec*scale[kax]
                        v11 = (Vals[0][1][ks]-Vals[2][1][ks]*l2)*vec*scale[kax]
                        M1.append(Vals[0][1][ks]+Vals[2][1][ks]*l2)
                        M1.append(Vals[0][1][ks]-Vals[2][1][ks]*l2)
                        ax.plot([x[0]+v00[0],x[0]+v01[0],x[1]+v11[0],x[1]+v10[0],x[0]+v00[0]],
                                [y[0]+v00[1],y[0]+v01[1],y[1]+v11[1],y[1]+v10[1],y[0]+v00[1]],'k',linewidth=0.5)
                        ax.annotate('%1.0f'%(Vals[kax][0][ks]),
                                    xy=(0.5*(x[0]+x[1])+v00[0],0.5*(y[0]+y[1])+v00[1]),
                                    ha=ha0,va='center',size=8)
                        ax.annotate('%1.0f'%(Vals[kax][1][ks]),
                                    xy=(0.5*(x[0]+x[1])+v01[0],0.5*(y[0]+y[1])+v01[1]),
                                    ha=ha1,va='center',size=8)
                    else:
                        v0 = Vals[kax][0][ks]*vec*scale[kax]
                        v1 = Vals[kax][1][ks]*vec*scale[kax]
                        ax.plot([x[0]+v0[0],x[0]+v1[0],x[1]+v1[0],x[1]+v0[0],x[0]+v0[0]],
                                [y[0]+v0[1],y[0]+v1[1],y[1]+v1[1],y[1]+v0[1],y[0]+v0[1]],'k',linewidth=0.5)
                        ax.annotate('%1.0f'%(Vals[kax][0][ks]),
                                    xy=(0.5*(x[0]+x[1])+v0[0],0.5*(y[0]+y[1])+v0[1]),
                                    ha=ha0,va='center',size=8)
                        ax.annotate('%1.0f'%(Vals[kax][1][ks]),
                                    xy=(0.5*(x[0]+x[1])+v1[0],0.5*(y[0]+y[1])+v1[1]),
                                    ha=ha1,va='center',size=8)
                    ax.plot([-33.07,33.07],[44,44],'k:',linewidth=0.5)
                    ax.plot([-33.07,33.07],[36.15,36.15],'k:',linewidth=0.5)
                    ax.plot([-33.07,33.07],[28.47,28.47],'k:',linewidth=0.5)
                    ax.plot([-33.07,33.07],[16.8,16.8],'k:',linewidth=0.5)
    ##    for kk in range(1,3):
    ##        minmax[kk][0] = min(minmax[kk][0],min(Vals[kk]))
    ##        minmax[kk][1] = max(minmax[kk][1],max(Vals[kk]))
    ##    minmax[0][0] = min(minmax[0][0],min(M))
    ##    minmax[0][1] = max(minmax[0][1],max(M))
    ##
        AX[1].annotate('Moment usuel\n[kNm/m]',
                       xy=(0.5,0.3),xycoords='axes fraction',ha='center',va='center')
        AX[2].annotate('Effort normal\n[kN/m]',
                       xy=(0.5,0.3),xycoords='axes fraction',ha='center',va='center')
        AX[3].annotate('Effort tranchant\n[kN/m]',
                       xy=(0.5,0.3),xycoords='axes fraction',ha='center',va='center')

        AX[0].set_title(u'Enveloppe des efforts, coupe à %1.0f° p.r. à\nl\'axe des tunnels'%(180/math.pi*(zvect[kz]+math.pi/2)))

        for seg in outline_segments:
            AX[0].plot(seg[0],seg[1],'k')
        n = plane.GetNormal()
        R = 33.07+5
        AX[0].plot([R*n[2]+n[0]*5,R*n[2],-R*n[2],-R*n[2]+n[0]*5],
                   [-R*n[0]+n[2]*5,-R*n[0],R*n[0],R*n[0]+n[2]*5],'b')

        for ax in AX[1:]:
    ##        ax.set_xlim([-6.4,0.4])
    ##        ax.set_ylim([-2.2,8.2])
            ax.set_ylabel('Altitude [NGF]')
    ##        ax.set_aspect('equal')
        ##    ax.set_title('Coupe '+titstr[k],size=20)
        AX[0].set_aspect('equal')

            
        fig.tight_layout()
        fig.savefig(prob+'_cut_Envelope_%i'%(kz))
        plt.close(fig)

    



        

            




