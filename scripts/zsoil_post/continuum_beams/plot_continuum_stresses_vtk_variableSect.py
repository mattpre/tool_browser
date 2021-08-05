# -*- coding: cp1252 -*-
# @description Plotting stresses on volumic cross sections.
# @input zsoil results
# @output png
# @author Matthias Preisig
# @date 2018/10/24

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import vtk
import os

from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools


loc_syst = np.array([[1,0,0],
                     [0,1,0],
                     [0,0,1]])

pathname = '../profil001'
pblist = ['M1161_P001_v2']

tvect = [4,5,7]

stepdict = {4:u'état existant',
            5:u'excavation du tunnel Ouest',
            7:u'excavation du tunnel Est'}

tx = []
for prob in pblist:
    try:
        f = open('pv/'+prob+'.pvd')
        for line in f:
            if 'vol.vtu' in line:
                t = float(line[line.find('timestep=')+9:].split('"')[1])
                if t in tvect:
                    tx.append(t)
        f.close()
        print('pv/'+prob+'.pvd loaded')
    except:
        res = zr(pathname,prob)
        res.read_rcf()
        res.read_his()
        tsteps = []
        for kt,step in enumerate(res.steps):
            if step.conv_status in [-1]:
                if step.time in tvect:
                    tsteps.append(kt)
                    tx.append(step.time)
        res.out_steps = tsteps
        res.read_dat()
        res.read_s00()
        for lab in res.ele_group_labels:
            if lab=='VOLUMICS':
                res.read_s01()
        ptemp = os.getcwd()
        os.chdir('pv')
        vtktools.write_vtu(res,vol=True,verbose=False)
        os.chdir(ptemp)

    # axe neutre tunnel LO:
    axe = []
    f = open('axe_LO_P001.txt')
    for line in f:
        v = line.split()
        axe.append((float(v[0]),float(v[1])))
    f.close()

    matlist = [6]
    sc = 0.0001
        
    for kt in range(len(tx)):
        t = tx[kt]

        tstr = '_'+str(int(t)).rjust(3,'0')+'_'+str(int((t-int(t))*100)).rjust(2,'0')
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName('pv/'+prob+tstr+'_vol.vtu')
        reader.Update()
        mesh = reader.GetOutput()
        
        MZ = []
        xN = []
        PT = []
        N = []
        h = []

        fig = plt.figure(figsize=(20,14))
        bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
        fig.text(0.01,0.01,prob,size=8)
        ax = fig.add_subplot(111)

        mat = mesh.GetCellData().GetArray('mat')
        patches = []
        cvect = []
        for kc in range(mesh.GetNumberOfCells()):
            cell = mesh.GetCell(kc)
            if mat.GetTuple1(kc) in matlist:
                pts = cell.GetPoints()
                xx = [pts.GetPoint(kp)[0] for kp in range(cell.GetNumberOfPoints())]
                yy = [pts.GetPoint(kp)[1] for kp in range(cell.GetNumberOfPoints())]
                patches.append(Polygon(np.array([xx,yy]).T))
                cvect.append(mat.GetTuple1(kc))
        pc = PatchCollection(patches,edgecolors='none',alpha=0.24,cmap=plt.cm.get_cmap('Set1'))
        pc.set_array(np.array(cvect))
        ax.add_collection(pc)
        
        for ka in range(len(axe)-1):
            patches = []
            cvect = []
            
            orig = np.array([0.5*(axe[ka+1][0]+axe[ka][0]),
                             0.5*(axe[ka+1][1]+axe[ka][1]),0])
            vect = np.array([(axe[ka+1][0]-axe[ka][0]),
                             (axe[ka+1][1]-axe[ka][1]),0])
            normal = vect/np.linalg.norm(vect)
            tang = np.array([-normal[1],normal[0],0])

            plane = vtk.vtkPlane()
            plane.SetOrigin(orig)
            plane.SetNormal(normal)
            
            cval,crd,co = vtktools.get_section_vol(mesh,plane,loc_syst,celldata=True,
                                                      array='STRESSES',matlist=matlist)
            cd = co.GetCellData()
            carr = cd.GetArray('STRESSES')
            pd = co.GetPointData()
            parr = pd.GetArray('DISP_TRA')

            r = []
            for ke in range(co.GetNumberOfCells()):
                cell = co.GetCell(ke)
                for kk in range(cell.GetNumberOfPoints()):
                    pt = co.GetPoint(cell.GetPointId(kk))
                    v = np.dot(np.array(pt)-orig,tang)
                    if abs(v)<2:
                        r.append(v)

            if len(r):
                rmin = min(r)
                rmax = max(r)
                h.append(rmax-rmin)
                pt0 = orig + 0.5*(rmax+rmin)*tang

                Mz = 0
                Nx = 0
                slist = []

                for ke in range(co.GetNumberOfCells()):
                    cell = co.GetCell(ke)
                    p0 = np.array(co.GetPoint(cell.GetPointId(0)))
                    p1 = np.array(co.GetPoint(cell.GetPointId(1)))
                    ptm = 0.5*(p0+p1)
                        
                    v = np.dot(np.array(ptm)-orig,tang)
                    if abs(v)<2:
                        stress = carr.GetTuple(ke)

                        length = cell.GetLength2()**0.5
                        xx = np.dot(ptm-pt0,tang)
                        yy = np.dot(ptm-pt0,normal)

                        strt0 = np.array([[stress[0],stress[2],0],
                                         [stress[2],stress[1],0],
                                         [0,0,stress[3]]])
                        cs = normal[0]
                        ss = normal[1]
                        rot = np.array([[cs,ss,0],
                                        [-ss,cs,0],
                                        [0,0,1]])
                        strt2 = np.matmul(rot,strt0)
                        strt = np.matmul(strt2,rot.T)
                        
                        Mz += length*xx*strt[0,0]
                        Nx += length*strt[0,0]

                        if ka%4==0:
                            p01 = p0 + normal*sc*strt[0,0]
                            p11 = p1 + normal*sc*strt[0,0]
                            xx = [p0[0],p01[0],p11[0],p1[0]]
                            yy = [p0[1],p01[1],p11[1],p1[1]]
                            patches.append(Polygon(np.array([xx,yy]).T))
                            if strt[0,0]<0:
                                cvect.append('b')
                            else:
                                cvect.append('r')
                            slist.append(strt[0,0])

                MZ.append(Mz)
                xN.append(Mz/Nx)
                PT.append(pt0)
                N.append(tang)
                if ka%4==0:
                    ax.annotate('$\sigma_{max}=$%1.2f MPa\n$\sigma_{min}=$%1.2f MPa'%(max(slist)*1e-3,min(slist)*1e-3),
                                xy=(pt0[0]-1.5*tang[0],pt0[1]-tang[1]),
                                ha='center',va='center',size=14)

            pc = PatchCollection(patches,edgecolors='none',facecolors=cvect,cmap=plt.cm.get_cmap('jet'),antialiased=True)
##            pc.set_array(np.array(cvect))
            ax.add_collection(pc)

        scm = 1./max(max(MZ),-min(MZ))
        ax.plot([PT[kp][0]+N[kp][0]*scm*MZ[kp] for kp in range(len(PT))],
                [PT[kp][1]+N[kp][1]*scm*MZ[kp] for kp in range(len(PT))])
        ax.plot([PT[kp][0] for kp in range(len(PT))],
                [PT[kp][1] for kp in range(len(PT))],'k--')
        ax.plot([PT[kp][0]+N[kp][0]*1./6*h[kp] for kp in range(len(PT))],
                [PT[kp][1]+N[kp][1]*1./6*h[kp] for kp in range(len(PT))],'k-.')
        ax.plot([PT[kp][0]-N[kp][0]*1./6*h[kp] for kp in range(len(PT))],
                [PT[kp][1]-N[kp][1]*1./6*h[kp] for kp in range(len(PT))],'k-.')
        ax.plot([PT[kp][0]+N[kp][0]*xN[kp] for kp in range(len(PT))],
                [PT[kp][1]+N[kp][1]*xN[kp] for kp in range(len(PT))],'r')
        fig.text(0.5,0.35,
                 u'Etape T=%i: %s'%(t,stepdict[t]),
                 ha='center',va='top',size=20)
        fig.text(0.5,0.32,
                 u'Moments dans la maçonnerie\n$M_{min}=$%1.1f kNm\n$M_{max}=$%1.1f kNm'%(min(MZ),max(MZ)),
                 ha='center',va='top',size=16)
                
##        ax.set_xlim([-5.9,5.9])
##        ax.set_ylim([-0.1,8])
        ax.axis('off')
        ax.set_aspect('equal')
##        ax.annotate(titstr[kax],xy=(0,2),ha='center',bbox=dict(pad=5,fc='w',ec='k'))
##            ax.set_title(titstr[kax])
##                xl = ax.get_xlim()
##                ax.plot(xl,[36.15,36.15],'k:')
##                ax.plot(xl,[28.47,28.47],'k:')
##                ax.plot(xl,[22.3,22.3],'k-.')
##                ax.plot(xl,[16.8,16.8],'k:')
##                ax.set_xlim(xl)
        
                

                
        fig.tight_layout()
        fig.savefig(prob+'_stresses_T%i'%(t))
        plt.close(fig)

    



        

            




