# -*- coding: cp1252 -*-
# @description settlement map and settlement profiles along sections. The sections are outlined on the map.
# @input zsoil results
# @output png
# @author Matthias Preisig
# @project M1147 SDP
# @date 2018/10/23

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import vtk
import os

from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools

batiments = [[(-6.26407e+01,2.90326e+01),
              (-3.62187e+01,2.84394e+01),
              (-3.58138e+01,5.56651e+01),
              (-6.24389e+01,5.60611e+01)],
             [(3.22243e+01,3.94354e+01),
              (3.47695e+01,4.90088e+01),
              (6.35670e+01,4.65781e+01),
              (6.28235e+01,3.69780e+01)],
             [(6.69553e+01,4.55266e+01),
              (6.74390e+01,5.11259e+01),
              (8.58266e+01,4.88501e+01),
              (8.55461e+01,4.34480e+01)],
             [(9.26095e+01,3.09931e+01),
              (9.40950e+01,4.50174e+01),
              (1.39951e+02,3.96427e+01),
              (1.38481e+02,2.48847e+01)]]
batlabs = ['4 (R+8)','5 (R+1)','6 (R+1)','7-8 (R+2)']


pathname = '..'
pblist = ['M1147_3D_SDP_coarse_v3_5(creep)_vo_18_0',
          'M1147_3D_SDP_coarse_v3_5(creep)_vo_315576000_0']
tx = []
X = [[] for prob in pblist]
UUY = [[] for prob in pblist]
for kpb,prob in enumerate(pblist):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName('../pv/'+prob+'.vtu')
    reader.Update()
    mesh = reader.GetOutput()
    
    plane = vtk.vtkPlane()
    plane.SetOrigin([0,37.5,0])
    plane.SetNormal([0,1,0])
    loc_syst = np.array([(1,0,0),(0,1,0)])
    val_y,crd,output = vtktools.get_section_vol(mesh,plane,loc_syst,
                                                array='disp',component=1)

    
    plane_xz = vtk.vtkPlane()
    plane_xz.SetOrigin([0,37.4,0])
    plane_xz.SetNormal([0,1,0])
    loc_syst_xz = np.array([(1,0,0),(0,0,-1)])

    val_y,crd,output = vtktools.get_section_vol(mesh,plane_xz,loc_syst_xz,
                                                array='disp',component=1)
##    w = vtk.vtkPolyDataWriter()
##    w.SetFileName('cut_%i.vtk'%(k))
##    w.SetInputData(output)
##    w.Write()

    fig_xz = plt.figure(figsize=(14,10))
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    fig_xz.text(0.01,0.01,prob,size=8)

    ax_xz = fig_xz.add_subplot(111)
    val1 = [max(min(0,val_y[kv]),-0.03001)*1e3 for kv in range(len(val_y))]

    CS = vtktools.contourf(val1,crd,output,loc_syst_xz,[0,37.4,0],ax_xz,levels=11)
    
    cb = fig_xz.colorbar(CS)
    cb.set_label('Tassement [mm]',size=12)
    cb.ax.tick_params(labelsize=12)

    # plot contours bâtiments:
    oo = []
    for kb in range(len(batiments)):
        ax_xz.plot([batiments[kb][kk][0] for kk in [0,1,2,3,0]],
                   [-batiments[kb][kk][1] for kk in [0,1,2,3,0]],'k')
        oo.append([np.mean([bb[0] for bb in batiments[kb]]),0,
                   np.mean([bb[1] for bb in batiments[kb]])])

    ##############################################################################
    # coupes:
    ##############################################################################
    titstr = ['A - A',
              'B - B']
    orig = [[1.28583e+02,0,-4.52748e+01],
            [1.28583e+02,0,-4.52748e+01]]
    normal = [[48.95,0,-14.147],
              [-14.147,0,-48.95]]
    BOUNDS = [[-80,80],
              [-10,100]]
    # add coupes par les centres des bâtiments:
    for kb in range(len(batiments)):
        titstr.append(batlabs[kb])
        orig.append(oo[kb])
        if kb==0:
            normal.append([1,0,1])
        else:
            normal.append([1,0,0])
        BOUNDS.append([-80,40])

    for k in range(len(orig)):
        plane = vtk.vtkPlane()
        plane.SetOrigin(orig[k])
        plane.SetNormal(normal[k])
        base = np.array([np.cross(normal[k]/np.linalg.norm(normal[k]),(0,1,0)),(0,1,0)])

        cutEdges = vtk.vtkCutter()
        cutEdges.SetInputData(output)
        cutEdges.SetCutFunction(plane)
        cutEdges.SetGenerateTriangles(0)
        cutEdges.Update()
    
        cuvette = cutEdges.GetOutput()
        points = cuvette.GetPoints()
        pdata = cuvette.GetPointData()
        UY = pdata.GetArray('disp')

        cuv_data = []
        pts = []
        for kpt in range(points.GetNumberOfPoints()):
            pt0 = points.GetPoint(kpt)
##            pts[1].append(orig[k][2]-pt0[2])
            pt1 = vtktools.project_on_plane(base,orig[k],pt0)
            if pt1[0]>BOUNDS[k][0] and pt1[0]<BOUNDS[k][1]:
                cuv_data.append([pt1[0],UY.GetTuple(kpt)[1]])
                pts.append([pt0[0],-pt0[2]])
        cuv_data.sort(key = lambda v:v[0])
        X[kpb].append([v[0] for v in cuv_data])
        UUY[kpb].append([v[1]*1e3 for v in cuv_data])


        pts.sort(key = lambda v:v[0])
        ax_xz.plot([pts[0][0],pts[-1][0]],
                   [pts[0][1],pts[-1][1]],'k')



##                    ax.annotate('Min.: %1.0f mm'%(min([val1[kk] for kk in range(len(crd[0])) if crd[0][kk]>40])),
##                                xy=(65,530),bbox=dict(fc='w',ec='k'))
        ang = np.angle((pts[-1][0]-pts[0][0])+(pts[-1][1]-pts[0][1])*1j)/np.pi*180
##        ax_xz.annotate(titstr[k],xy=pts[-5],bbox=dict(fc='w',ec='k'),rotation=ang)
        ax_xz.annotate(titstr[k].split()[0],xy=pts[0],rotation=ang)
        ax_xz.annotate(titstr[k].split()[0],xy=pts[-1],rotation=ang)

    ax_xz.set_xlim([-110,230])
    ax_xz.set_ylim([-100,180])
    ax_xz.set_xlabel(u'Est - Ouest',size=12)
    ax_xz.set_ylabel(u'Nord - Sud',size=12)
    ax_xz.set_aspect('equal')

        
##            ax.set_title('Coupe '+titstr[k],size=20)
        
    fig_xz.tight_layout()
    fig_xz.savefig(prob+'_cuvetteLines')
    plt.close(fig_xz)


labstr = ['Fin des travaux',
          'Long terme']
crit = [[15,8],
        [3,7],
        [5,7.5],
        [5,7.5],
        [5,7.5],
        [5,7.5]]
cpos = [[20,-17],
        [3,-19],
        [2.5,-2],
        [2.5,-10],
        [2.5,-10],
        [2.5,-10]]
for k in range(len(orig)):
    fig = plt.figure(figsize=(10,7))
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    fig.text(0.01,0.01,prob,size=8)

    ax = fig.add_subplot(111)

    for kpb,prob in enumerate(pblist):
        ax.plot(X[kpb][k],UUY[kpb][k],label=labstr[kpb])

##    ax.set_xlim(BOUNDS[k])
##        ax.set_ylim([-2.2,8.2])
    ax.set_ylabel('Tassement [mm]')
    ax.legend()

    ax.plot([cpos[k][0]-crit[k][0],
             cpos[k][0],
             cpos[k][0],
             cpos[k][0]-crit[k][0]],
            [cpos[k][1],
             cpos[k][1],
             cpos[k][1]+crit[k][1],
             cpos[k][1]],'k',linewidth=2)
    
    ax.grid('on')
    ax.set_title('Coupe '+titstr[k],size=20)

        
    fig.tight_layout()
    fig.savefig(prob+'_cuvette_'+titstr[k].split()[0])
    plt.close(fig)


            




