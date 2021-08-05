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
import re

from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools

#           nom coupe   niveau   origine          dir. coupe   range    pos point
cuvettes = [['Coupe 1',496,(0.5,496,0),(0,0,1),(-40,40),[]],
            ['Coupe 2',496,(20,496,0),(0,0,1),(-40,40),[]]]

niveaux = list(set([v[1] for v in cuvettes]))

pathname = '..'
pblist = ['M1161_caverne_v1_2']
tvect = [2,11,12,24,35,41,52]
stepdict = {2:u'excavation galerie en calotte',
            11:u'poutre et dalle',
            12:u'approfondissement stross',
            21:u'jetting tunnel, voûte et front',
            22:u'démolition radier intersection',
            23:u'excavation et réalisation radier tunnel',
            24:u'ouverture galerie côté montagne',
            40:'nn',
            41:u'ouverture galerie côté lac'}
for kk in range(12):
    stepdict[25+kk] = u'front d\'excavation montagne à %1.0fm'%(kk+1)
    stepdict[42+kk] = u'front d\'excavation lac à %1.0fm'%(kk+1)
t0 = 2
tx = []
X = [[] for t in tvect]
UUY = [[] for t in tvect]
UUH = [[] for t in tvect]
for kpb,prob in enumerate(pblist):
    meshes = []
    try:
        for t in tvect:
            tstr = vtktools.get_tstr(t)
            open('pv/'+prob+'_'+tstr+'_vol.vtu').close()
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
##                elif abs(step.time-t0)<1e-3:
##                    tsteps.append(kt)
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
                if step.time in tvect:
                    res.out_steps.append(kt)
                elif step.time==t0:
                    step0 = res.steps[kt]
        vtktools.write_vtu(res,vol=True,verbose=False,refstep=step0)
        os.chdir(ptemp)

    for kt,t in enumerate(tvect):
        reader = vtk.vtkXMLUnstructuredGridReader()
        tstr = vtktools.get_tstr(t)
        print('opening file %s'%('pv/'+prob+'_'+tstr+'_vol.vtu'))
        reader.SetFileName('pv/'+prob+'_'+tstr+'_vol.vtu')
        reader.Update()
        mesh = reader.GetOutput()

    ##    w = vtk.vtkPolyDataWriter()
    ##    w.SetFileName('cut_%i.vtk'%(1))
    ##    w.SetInputData(output)
    ##    w.Write()

        planes = []
        for kniv,niv in enumerate(niveaux):
            plane = vtk.vtkPlane()
            plane.SetOrigin([0,niv,0])
            plane.SetNormal([0,1,0])
            loc_syst = np.array([(1,0,0),(0,1,0)])
            val_y,crd,output = vtktools.get_section_vol(mesh,plane,loc_syst,
                                                        array='DISP_TRA',component=1)
            planes.append(output)
            
        ##############################################################################
        # coupes:
        ##############################################################################

        for kc,coupe in enumerate(cuvettes):
            plane = vtk.vtkPlane()
            plane.SetOrigin(coupe[2])
            normal = np.array([coupe[3][2],0,-coupe[3][0]])
            normal = normal/np.linalg.norm(normal)
            plane.SetNormal(normal)
            base = np.array([np.cross(normal,(0,1,0)),(0,1,0)])

            ind = niveaux.index(coupe[1])

            cutEdges = vtk.vtkCutter()
            cutEdges.SetInputData(planes[ind])
            cutEdges.SetCutFunction(plane)
            cutEdges.SetGenerateTriangles(0)
            cutEdges.Update()
        
            cuvette = cutEdges.GetOutput()
            points = cuvette.GetPoints()
            pdata = cuvette.GetPointData()
            UY = pdata.GetArray('DISP_TRA')
            
    ##        w = vtk.vtkPolyDataWriter()
    ##        w.SetFileName('cuvette_%i.vtk'%(kc))
    ##        w.SetInputData(cuvette)
    ##        w.Write()

            cuv_data = []
            pts = []
            for kpt in range(points.GetNumberOfPoints()):
                pt0 = points.GetPoint(kpt)
                pt1 = vtktools.project_on_plane(base,coupe[2],pt0)
                if pt1[0]>coupe[4][0] and pt1[0]<coupe[4][1]:
                    uh = np.dot(base[0],np.array([UY.GetTuple(kpt)[0],0,UY.GetTuple(kpt)[2]]))
                    cuv_data.append([pt1[0],UY.GetTuple(kpt)[1],uh])
                    pts.append(([pt0[0],-pt0[2]],pt1[0]))
            cuv_data.sort(key = lambda v:v[0])
            pts.sort(key = lambda v:v[1])
            pts = [v[0] for v in pts]
            X[kt].append([v[0] for v in cuv_data])
            UUY[kt].append([v[1]*1e3 for v in cuv_data])
            UUH[kt].append([v[2]*1e3 for v in cuv_data])


for k in range(len(cuvettes)):
    fig = plt.figure(figsize=(12,10))
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    fig.text(0.01,0.01,prob,size=8)

    ax = fig.add_subplot(111)

    for kt,t in enumerate(tvect):
        ax.plot(X[kt][k],UUY[kt][k],label=stepdict[t])

    ax.set_xlim([-17,17])
    ax.set_ylim([ax.get_ylim()[0],0])
    ax.set_ylabel('Tassement en surface [mm]')
    ax.set_xlabel('Lac'+' '*20+'<->'+' '*20+'Montagne')
    ax.legend()
    
    ax.grid('on')
    ax.set_title(cuvettes[k][0],size=20)

        
    fig.tight_layout()
    fig.savefig(prob+'_cuv_'+re.sub(' ','_',cuvettes[k][0]))
    plt.close(fig)



