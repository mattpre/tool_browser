# -*- coding: cp1252 -*-
# @description Plotting M, N, T-diagrams from integrated volumic elements along cross sections.
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

pathname = '..'
pblist = ['M1135_noISS_v2_3']
pblist = ['M1135_Var1_v3_2(pseudoEnka)']
##pblist = ['M1135_Var2_v3_1']

tvect = [7,8,10,14,15,17]
t0vect = [6,6,6,12,12,12]
t0 = 6

matdict = {1:u'rocher',
           2:u'béton ancien',
           3:u'béton neuf (coffré)',
           4:u'béton projeté',
           5:u'étanchéité'}
stepdict = {7:u'clouage et renforcement voûte',
            8:u'action bloc local sur ancien béton',
            10:u'action uniforme sur ancien béton',
            13:u'revêtement posé',
            14:u'clous désactivés',
            15:u'action bloc local sur revêtement',
            17:u'action uniforme sur revêtement'}

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
        tx = tvect
        print('pv/'+prob+'.pvd loaded')
    except:
        res = zr(pathname,prob)
        res.read_rcf()
        res.read_his()
        tsteps = []
        step0 = False
        for kt,step in enumerate(res.steps):
            if step.conv_status in [-1]:
                if step.time in tvect:
                    tsteps.append(kt)
                    tx.append(step.time)
                elif step.time in t0vect:
                    tsteps.append(kt)
        res.out_steps = tsteps
        res.read_dat()
        res.read_s00()
        for lab in res.ele_group_labels:
            if lab=='VOLUMICS':
                res.read_s01()
        ptemp = os.getcwd()
        os.chdir('pv')
        res.out_steps = []
        for kt,step in enumerate(res.steps):
            if step.conv_status in [-1]:
                if step.time in [7,8,10]:
                    res.out_steps.append(kt)
                elif step.time==6:
                    step0 = res.steps[kt]
        vtktools.write_vtu(res,vol=True,verbose=False,refstep=step0)
        res.out_steps = []
        for kt,step in enumerate(res.steps):
            if step.conv_status in [-1]:
                if step.time in [14,15,17]:
                    res.out_steps.append(kt)
                elif step.time==12:
                    step0 = res.steps[kt]
        vtktools.write_vtu(res,vol=True,verbose=False,refstep=step0)
        os.chdir(ptemp)

        
    for kt in range(len(tx)):
        t = tx[kt]
        t0 = t0vect[kt]

        if t<12:
            matlist = [2,4]
            a0 = math.pi*0.64
        else:
            matlist = [3]
            a0 = math.pi*0.57
        nA = 401
        avect = [2*a0/(nA-1)*k-a0 for k in range(nA)]

        tstr = '_'+str(int(t)).rjust(3,'0')+'_'+str(int((t-int(t))*100)).rjust(2,'0')
        tstr += '-'+str(int(t0)).rjust(3,'0')+'_'+str(int((t0-int(t0))*100)).rjust(2,'0')
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName('pv/'+prob+tstr+'_vol.vtu')
        reader.Update()
        mesh = reader.GetOutput()

        MZ = []
        NX = []
        QZ = []
        UN = []
        U = []
        PT = []
        N = []
        
        orig = np.array([0,2,0])
        for ka in range(len(avect)):
            normal = np.array([math.cos(avect[ka]),-math.sin(avect[ka]),0])
            tang = np.array([math.sin(avect[ka]),math.cos(avect[ka]),0])

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
                    if v>0:
                        r.append(v)

            if len(r):
                rmin = min(r)
                rmax = max(r)
                pt0 = orig + 0.5*(rmax+rmin)*tang

                L = 0
                Mz = 0
                Nx = 0
                Qz = 0
                ddist = 1e10
                un = 0
                
                for ke in range(co.GetNumberOfCells()):
                    cell = co.GetCell(ke)
                    ptm = np.array(co.GetPoint(cell.GetPointId(0)))
                    ptm += np.array(co.GetPoint(cell.GetPointId(1)))
                    ptm *= 0.5
                        
                    v = np.dot(np.array(ptm)-orig,tang)
                    if v>0:
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

                        dist = np.linalg.norm(pt0-co.GetPoint(cell.GetPointId(0)))
                        if dist<ddist:
                            disp = np.array(parr.GetTuple(cell.GetPointId(0)))
                            rot2 = np.array([[cs,-ss],
                                            [ss,cs]])
                            drt = np.dot(rot2,disp)
                            ddist = dist
                        dist = np.linalg.norm(pt0-co.GetPoint(cell.GetPointId(1)))
                        if dist<ddist:
                            disp = np.array(parr.GetTuple(cell.GetPointId(1)))
                            rot2 = np.array([[cs,-ss],
                                            [ss,cs]])
                            drt = np.dot(rot2,disp)
                            ddist = dist
                        
                        L += length
                        Mz += length*xx*strt[0,0]
                        Nx += length*strt[0,0]
                        Qz += length*strt[0,1]
                        un = drt[1]
    ##            x0 = [Ax[kk]/L[kk] for kk in range(len(L))]
    ##            z0 = [Az[kk]/L[kk] for kk in range(len(L))]
    ##            ang = [math.atan(z0[kk]/x0[kk]) for kk in range(len(Area))]
    ##            Mz = [Mz0[kk]-(Ny[kk]*z0[kk]) for kk in range(len(Area))]

                MZ.append(Mz)
                NX.append(Nx)
                QZ.append(Qz)
                PT.append(pt0)
                N.append(tang)
                UN.append(un*1000)
                U.append(disp*1000)

        fig = plt.figure(figsize=(14,11))
        bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
        fig.text(0.01,0.01,prob,size=8)
        matstr = str()
        for km in matlist:
            matstr += ', %s'%(matdict[km])
        fig.text(0.5,0.98,
                 u'Etape T=%i: %s\nEfforts dans les matériaux: %s'%(t,stepdict[t],matstr[2:]),
                 ha='center',va='top',size=16)

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


        scn = 1./max(max(NX),-min(NX))
        scm = 1./max(max(MZ),-min(MZ))
        scq = 1./max(max(QZ),-min(QZ))
        scu = 1./max(max(UN),-min(UN))
        AX = []
        ax = fig.add_subplot(221)
        ax.plot([PT[kp][0]+N[kp][0]*scn*NX[kp] for kp in range(len(PT))],
                [PT[kp][1]+N[kp][1]*scn*NX[kp] for kp in range(len(PT))])
        kmax = NX.index(max(NX))
        kmin = NX.index(min(NX))
        ax.plot(PT[kmax][0]+N[kmax][0]*scn*NX[kmax],
                PT[kmax][1]+N[kmax][1]*scn*NX[kmax],'g.')
        ax.plot(PT[kmin][0]+N[kmin][0]*scn*NX[kmin],
                PT[kmin][1]+N[kmin][1]*scn*NX[kmin],'r.')
        ax.plot([PT[kp][0] for kp in range(len(PT))],
                [PT[kp][1] for kp in range(len(PT))],'k')
        pc = PatchCollection(patches,edgecolors='none',alpha=0.2,cmap=plt.cm.get_cmap('Set1'))
        pc.set_array(np.array(cvect))
        ax.add_collection(pc)
        AX.append(ax)
        ax = fig.add_subplot(222)
        ax.plot([PT[kp][0]+N[kp][0]*scm*MZ[kp] for kp in range(len(PT))],
                [PT[kp][1]+N[kp][1]*scm*MZ[kp] for kp in range(len(PT))])
        kmax = MZ.index(max(MZ))
        kmin = MZ.index(min(MZ))
        ax.plot(PT[kmax][0]+N[kmax][0]*scm*MZ[kmax],
                PT[kmax][1]+N[kmax][1]*scm*MZ[kmax],'g.')
        ax.plot(PT[kmin][0]+N[kmin][0]*scm*MZ[kmin],
                PT[kmin][1]+N[kmin][1]*scm*MZ[kmin],'r.')
        ax.plot([PT[kp][0] for kp in range(len(PT))],
                [PT[kp][1] for kp in range(len(PT))],'k')
        pc = PatchCollection(patches,edgecolors='none',alpha=0.2,cmap=plt.cm.get_cmap('Set1'))
        pc.set_array(np.array(cvect))
        ax.add_collection(pc)
        AX.append(ax)
        ax = fig.add_subplot(223)
        ax.plot([PT[kp][0]+N[kp][0]*scq*QZ[kp] for kp in range(len(PT))],
                [PT[kp][1]+N[kp][1]*scq*QZ[kp] for kp in range(len(PT))])
        if max(QZ)>-min(QZ):
            kmax = QZ.index(max(QZ))
        else:
            kmax = QZ.index(min(QZ))
        ax.plot(PT[kmax][0]+N[kmax][0]*scq*QZ[kmax],
                PT[kmax][1]+N[kmax][1]*scq*QZ[kmax],'r.')
        ax.plot([PT[kp][0] for kp in range(len(PT))],
                [PT[kp][1] for kp in range(len(PT))],'k')
        pc = PatchCollection(patches,edgecolors='none',alpha=0.2,cmap=plt.cm.get_cmap('Set1'))
        pc.set_array(np.array(cvect))
        ax.add_collection(pc)
        AX.append(ax)
        ax = fig.add_subplot(224)
        ax.plot([PT[kp][0]+scu*U[kp][0] for kp in range(len(PT))],
                [PT[kp][1]+scu*U[kp][1] for kp in range(len(PT))])
        UM = [np.linalg.norm(u) for u in U]
        kmax = UM.index(max(UM))
        kmin = UM.index(min(UM))
        ax.plot(PT[kmax][0]+scu*U[kmax][0],
                PT[kmax][1]+scu*U[kmax][1],'g.')
        ax.plot(PT[kmin][0]+scu*U[kmin][0],
                PT[kmin][1]+scu*U[kmin][1],'r.')
        ax.plot([PT[kp][0] for kp in range(len(PT))],
                [PT[kp][1] for kp in range(len(PT))],'k')
        pc = PatchCollection(patches,edgecolors='none',alpha=0.2,cmap=plt.cm.get_cmap('Set1'))
        pc.set_array(np.array(cvect))
        ax.add_collection(pc)
        AX.append(ax)

        titstr = ['Effort normal\n$N_{min}=$%1.1f kN\n$N_{max}=$%1.1f kN'%(min(NX),max(NX)),
                  'Moment de flexion\n$M_{min}=$%1.1f kNm ($N=$%1.1f kN)\n$M_{max}=$%1.1f kNm ($N=$%1.1f kN)'%(min(MZ),NX[MZ.index(min(MZ))],
                                                                                                               max(MZ),NX[MZ.index(max(MZ))]),
                  'Effort tranchant\n$|Q|_{max}=$%1.1f kN'%(max(-min(QZ),max(QZ))),
                  u'Déplacement normal\n$u_{n,min}=$%1.2f mm\n$u_{n,max}=$%1.2f mm'%(min(UN),max(UN))]
        xlabstr = ['NX [kN]','MZ [kNm]','QZ [kN]']
        for kax,ax in enumerate(AX):
##            ax.grid(True)
            ax.set_xlim([-6.2,6.2])
            ax.set_ylim([-0.3,8.2])
##                ax.set_ylabel('Altitude [NGF]')
##                ax.set_xlabel(xlabstr[kax])
            ax.axis('off')
            ax.set_aspect('equal')
            ax.annotate(titstr[kax],xy=(0,2),ha='center',bbox=dict(pad=5,fc='w',ec='k'))
##            ax.set_title(titstr[kax])
##                xl = ax.get_xlim()
##                ax.plot(xl,[36.15,36.15],'k:')
##                ax.plot(xl,[28.47,28.47],'k:')
##                ax.plot(xl,[22.3,22.3],'k-.')
##                ax.plot(xl,[16.8,16.8],'k:')
##                ax.set_xlim(xl)
        
                

                
        fig.tight_layout()
        fig.savefig(prob+'_diags_T%i'%(t))
        plt.close(fig)

    



        

            




