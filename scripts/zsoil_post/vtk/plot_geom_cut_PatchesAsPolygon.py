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
import vtk.util.numpy_support as vnp
import os

from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools


loc_syst = np.array([[1,0,0],
                     [0,1,0],
                     [0,0,1]])

pathname = '..'
pblist = ['M1438_v2']

tvect = [1]

matdict = {1:u'Stahl',
           2:u'Beton'}

tx = []
for prob in pblist:
    try:
        f = open(pathname+'/pv/'+prob+'.pvd')
        for line in f:
            if 'vol.vtu' in line:
##            if '_vo_' in line:
                t = float(line[line.find('timestep=')+9:].split('"')[1])
                if t in tvect:
                    tx.append(t)
        f.close()
        print(pathname+'/pv/'+prob+'.pvd loaded')
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
            elif lab=='SHELLS':
                res.read_s02()
            elif lab=='CONTACTS':
                res.read_s07()
        ptemp = os.getcwd()
        os.chdir(pathname+'/pv')
        pvdFile = vtktools.create_pvd(res)
        vtktools.write_vtu(res,vol=True,verbose=False)
        vtktools.write_vtu(res,shells=True,verbose=False)
        vtktools.write_vtu(res,cnt=True,verbose=False)
        vtktools.save_pvd(pvdFile)
        os.chdir(ptemp)

    titstr = [u'Detail 1 Längsschnitt A Gurt unten',
              u'Detail 1 Schnitt B unten',
              u'Detail 1 Längsschnitt A Gurt oben',
              u'Detail 1 Schnitt C Gurt oben',
              u'Längsschnitt -5- Endquerträger durch Steg',
              u'Schnitt -2- Endquerträger',
              u'Schnitt Obergurt',
              u'Schnitt Druckstrebe',
              u'Schnitt Verankerung',
              u'Schnitt Untergurt']
    origs = [np.array([0,0,0.2]),
             np.array([5.9,0,0]),
             np.array([0,0,0.28]),
             np.array([4.1,0,0]),
             np.array([0,0,0.27]),
             np.array([5.9,0,0]),
             np.array([0,2.67,0]),
             np.array([3.366,2.67,0]),
             np.array([3.516,2.67,0]),
             np.array([0,0.93,0])]
    normals = [[0,0,1],
               [1,0,0],
               [0,0,1],
               [1,0,0],
               [0,0,1],
               [1,0,0],
               [0,1,0],
               [0.55455243,0.83214879,0],
               [-0.55685396,0.83061042,0],
               [0,1,0]]
    bases = [[np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])],
             [np.array([0,0,-1]),np.array([0,1,0]),np.array([1,0,0])],
             [np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])],
             [np.array([0,0,-1]),np.array([0,1,0]),np.array([1,0,0])],
             [np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])],
             [np.array([0,0,-1]),np.array([0,1,0]),np.array([1,0,0])],
             [np.array([1,0,0]),np.array([0,0,-1]),np.array([0,1,0])],
             [np.array([0.83214879,-0.55455243,0]),np.array([0,0,-1]),np.array([0.55455243,0.83214879,0])],
             [np.array([0.83061042,0.55685396,0]),np.array([0,0,-1]),np.array([-0.55685396,0.83061042,0])],
             [np.array([1,0,0]),np.array([0,0,-1]),np.array([0,1,0])]]
    bounds = [[[5.25,0.65],[7.6,1.4]],
              [[-0.52,0.7],[0.08,1.25]],
              [[2.7,1.7],[5,2.75]],
              [[-0.52,0.7],[0.08,2.75]],
              [[0.27,0.29],[7.25,2.75]],
              [[-0.52,0.78],[0.08,2.75]],
              [[2,-0.53],[5.0,0.08]],
              [[-0.2,-0.53],[4.0,0.08]],
              [[-4,-0.53],[0.2,0.08]],
              [[0.4,-0.53],[6,0.08]]]

    matlist = [1,2]

    t = 1

    tstr = '_'+str(int(t)).rjust(3,'0')+'_'+str(int((t-int(t))*100)).rjust(2,'0')
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(pathname+'/pv/'+prob+tstr+'_vol.vtu')
    reader.Update()
    volmesh = reader.GetOutput()

    cd = volmesh.GetCellData()
    carr = cd.GetArray('mat')

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(pathname+'/pv/'+prob+tstr+'_shell.vtu')
    reader.Update()
    shellmesh = reader.GetOutput()

    cd = shellmesh.GetCellData()
    sarr = cd.GetArray('mat')

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(pathname+'/pv/'+prob+tstr+'_cnt.vtu')
    reader.Update()
    cntmesh = reader.GetOutput()
    
    cmap=plt.colormaps['Set1']

    for kcoupe in range(len(titstr)):
        patches = []
        cvect = []

        lx = bounds[kcoupe][1][0]-bounds[kcoupe][0][0]
        ly = bounds[kcoupe][1][1]-bounds[kcoupe][0][1]        
        fig = plt.figure(figsize=(12/(ly+0.5)*lx,8))
        ax = fig.add_subplot(111)
        fig.text(0.01,0.01,prob,size=8)

        orig = origs[kcoupe]
        normal = normals[kcoupe]
        base = bases[kcoupe]

        plane = vtk.vtkPlane()
        plane.SetOrigin(orig)
        plane.SetNormal(normal)

        for mat in matlist:        
            cval,crd,co = vtktools.get_section_vol(volmesh,plane,base,celldata=True,
                                                   array='mat',matlist=[mat])

            cpd = vtk.vtkCleanPolyData()
            cpd.SetInputData(co)
            cpd.Update()
            
            if False:
                w = vtk.vtkPolyDataWriter()
                w.SetFileName('cut_%i.vtk'%(kcoupe))
                w.SetInputConnection(cpd.GetOutputPort())
                w.Write()
            fe = vtk.vtkFeatureEdges()
            fe.BoundaryEdgesOn()
            fe.NonManifoldEdgesOff()
            fe.ManifoldEdgesOff()
            fe.FeatureEdgesOff()
            fe.SetInputConnection(cpd.GetOutputPort())
            fe.Update()
            cpd = vtk.vtkCleanPolyData()
            cpd.SetInputConnection(fe.GetOutputPort())
            cpd.Update()
            co = cpd.GetOutput()
            
            if False:
                w = vtk.vtkPolyDataWriter()
                w.SetFileName('cut_ol_%i.vtk'%(kcoupe))
                w.SetInputData(co)
                w.Write()

    ##        print(co.GetBounds())
##            cd = co.GetCellData()
##            mat = cd.GetArray('mat')

            vertices = {}
            for kc in range(co.GetNumberOfCells()):
                cell = co.GetCell(kc)
                if cell.GetPointId(0) not in vertices:
                    vertices[cell.GetPointId(0)] = []
                if cell.GetPointId(1) not in vertices:
                    vertices[cell.GetPointId(1)] = []
                vertices[cell.GetPointId(0)].append(kc)
                vertices[cell.GetPointId(1)].append(kc)

##            vertexIds = [kp for kp in range(co.GetNumberOfPoints())]
            lines = [kc for kc in range(co.GetNumberOfCells())]
            while len(lines):
                kc = lines[0]
                polyIds = [co.GetCell(kc).GetPointId(0)]
                while len(polyIds)<1e5:
                    cell = co.GetCell(kc)
                    nextId = cell.GetPointId(1)
                    cellind = lines.index(kc)
                    lines.pop(cellind)
                    if nextId==polyIds[0]:
                        pts = [vtktools.project_on_plane(base,orig,co.GetPoint(kp)) for kp in polyIds]
                        patches.append(Polygon(pts))
                        cvect.append(mat)
                        polyIds = []
                        break
                    else:
                        polyIds.append(nextId)
                        ind = vertices[nextId].index(kc)
                        kc = vertices[nextId][(ind+1)%2]
##                    print(lines)
        print(len(cvect))

        pc = PatchCollection(patches,edgecolors='none',cmap=cmap,antialiased=False)
        pc.set_array(np.array(cvect))
        ax.add_collection(pc)
        
        # shells:
        segments = vtktools.get_section(shellmesh,plane,origin=orig,loc_syst=base,matlist=[5,8])
        for ks,seg in enumerate(segments):
            thickness = seg[2][3]
            p0 = seg[0]
            p1 = seg[1]
            n = np.array([p1[1]-p0[1],p0[0]-p1[0]])
            n /= np.linalg.norm(n)
            p00 = p0 - 0.5*n*thickness
            p01 = p0 + 0.5*n*thickness
            p10 = p1 - 0.5*n*thickness
            p11 = p1 + 0.5*n*thickness
            ax.add_patch(Polygon([p00,p01,p11,p10],ec='none',fc='b',alpha=0.5,antialiased=False))
##            ax.plot([seg[0][0],seg[1][0]],
##                    [seg[0][1],seg[1][1]],'k')

        # contacts:
        segments = vtktools.get_section(cntmesh,plane,origin=orig,loc_syst=base,matlist=[6])
        for ks,seg in enumerate(segments):
            p0 = seg[0]
            p1 = seg[1]
            ax.plot([seg[0][0],seg[1][0]],
                    [seg[0][1],seg[1][1]],'y',lw=5,alpha=0.8,solid_capstyle='butt')


##        cb = fig.colorbar(pc)
        ax.set_title(titstr[kcoupe])

        ax.set_xlim([bounds[kcoupe][0][0],bounds[kcoupe][1][0]])
        ax.set_ylim([bounds[kcoupe][0][1],bounds[kcoupe][1][1]])
##        ax.axis('off')
        ax.grid(axis='both')
        ax.set_aspect('equal')

        fig.tight_layout()
        fig.savefig(prob+'_geom_outline_%i.svg'%(kcoupe+1))
        plt.close(fig)

    



        

            




