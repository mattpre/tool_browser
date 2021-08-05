# -*- coding: cp1252 -*-
# @description Plotting stresses on volumic cross sections.
# @input zsoil results
# @output png
# @author Matthias Preisig
# @date 2018/10/24

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
import vtk
import os,re
from PIL import Image

from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools

lut = vtktools.get_lut(lut_type='mat',maxind=29)
clist = []
for k in range(1,21):
    col = lut.GetTableValue(k)
    clist.append(col)
cmap = colors.ListedColormap(clist)
norm = colors.BoundaryNorm(range(1,21),cmap.N)

loc_syst = np.array([[1,0,0],
                     [0,1,0],
                     [0,0,1]])

pathname = '../..'
pblist = ['M1161_P002_v1_1']

tvect = [0,1.2,2,3,4,4.99,5,5.7,5.99,5.99,6,6.7,7]

matdict = {1:u'Remblai',
           2:u'GL',
           3:u'Moraine de font',
           4:u'Molasse',
           6:u'Maçonnerie',
           7:u'Béton auge',
           8:u'Jetting',
           9:u'Béton remplissage',
           11:u'sout m2',
           12:u'sout m3',
           13:u'cnt',
           14:u'Hôtel Palace',
           15:u'Hôtel Palace'}
stepdict = {7:u'couage et renforcement voûte',
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
            if 'beam.vtu' in line:
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
            if lab=='TRUSSES':
                res.read_s03()
            if lab=='BEAMS':
                res.read_s04()
        ptemp = os.getcwd()
        os.chdir('pv')
        vtktools.write_vtu(res,vol=True,verbose=False)
        vtktools.write_vtu(res,beams=True,verbose=False)
        vtktools.write_vtu(res,trusses=True,verbose=False)
        os.chdir(ptemp)
        
    for kt in range(len(tx)):
        t = tx[kt]

        patches = []
        cvect = []
        lims = [[1e10,-1e10] for k in range(2)]

        tstr = vtktools.get_tstr(t)
        reader = vtk.vtkXMLUnstructuredGridReader()

        reader.SetFileName('pv/'+prob+'_'+tstr+'_vol.vtu')
        reader.Update()
        mesh = reader.GetOutput()

        mat = mesh.GetCellData().GetArray('mat')
        for kc in range(mesh.GetNumberOfCells()):
            cell = mesh.GetCell(kc)
##            if mat.GetTuple1(kc) in matlist:
            pts = cell.GetPoints()
            xx = [pts.GetPoint(kp)[0] for kp in range(cell.GetNumberOfPoints())]
            yy = [pts.GetPoint(kp)[1] for kp in range(cell.GetNumberOfPoints())]
            lims[0][0] = min(min(xx),lims[0][0])
            lims[0][1] = max(max(xx),lims[0][1])
            lims[1][0] = min(min(yy),lims[1][0])
            lims[1][1] = max(max(yy),lims[1][1])
            patches.append(Polygon(np.array([xx,yy]).T))
            cvect.append(mat.GetTuple1(kc))

        reader.SetFileName('pv/'+prob+'_'+tstr+'_beam.vtu')
        reader.Update()
        mesh = reader.GetOutput()

        lines = []
        cvect_beams = []
        mat = mesh.GetCellData().GetArray('mat')
        for kc in range(mesh.GetNumberOfCells()):
            cell = mesh.GetCell(kc)
##            if mat.GetTuple1(kc) in matlist:
            pts = cell.GetPoints()
            xx = [(pts.GetPoint(kp)[0],pts.GetPoint(kp)[1]) for kp in range(cell.GetNumberOfPoints())]
##            yy = [pts.GetPoint(kp)[1] for kp in range(cell.GetNumberOfPoints())]
            lines.append(xx)
            cvect_beams.append(mat.GetTuple1(kc))

        reader.SetFileName('pv/'+prob+'_'+tstr+'_truss.vtu')
        reader.Update()
        mesh = reader.GetOutput()

        mat = mesh.GetCellData().GetArray('mat')
        for kc in range(mesh.GetNumberOfCells()):
            cell = mesh.GetCell(kc)
##            if mat.GetTuple1(kc) in matlist:
            pts = cell.GetPoints()
            xx = [(pts.GetPoint(kp)[0],pts.GetPoint(kp)[1]) for kp in range(cell.GetNumberOfPoints())]
##            yy = [pts.GetPoint(kp)[1] for kp in range(cell.GetNumberOfPoints())]
            lines.append(xx)
            cvect_beams.append(mat.GetTuple1(kc))
        
        fig = plt.figure(figsize=(14,8))
        bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
        fig.text(0.01,0.01,prob,size=8)
        ax = fig.add_subplot(111)

        pc = PatchCollection(patches,edgecolors='k',cmap=cmap,norm=norm,
                             linewidth=0.5,antialiased=True)
        pc.set_array(np.array(cvect))
        ax.add_collection(pc)
        lc = LineCollection(lines,linewidths=3,cmap=cmap,norm=norm)
        lc.set_array(np.array(cvect_beams))
        ax.add_collection(lc)


##        fig.text(0.5,0.35,
##                 u'Etape T=%i: %s'%(t,stepdict[t]),
##                 ha='center',va='top',size=16)
                
        ax.set_xlim(lims[0])
        ax.set_ylim(lims[1])
        ax.axis('off')
        ax.set_aspect('equal')

                
        fig.tight_layout()
        fig.savefig(prob+'_phasage_%s'%(tstr))
        plt.close(fig)

        img = Image.open(prob+'_phasage_%s'%(tstr)+'.png')

        from zsoil_tools import postpro_lib as pl
        nval = lut.GetNumberOfColors()
        leg = pl.get_legend(lut,categories=matdict,hfrac=0.8,vpad=0.1,hwratio=4,
                            label='')
        
        leg = leg.resize((int(float(leg.size[0])/leg.size[1]*img.size[1]),img.size[1]))
        img1 = Image.new('RGB',(img.size[0]+leg.size[0],img.size[1]))
        img1.paste(img,(0,0))
        img1.paste(leg,(img.size[0],0))
        
        img1.save(prob+'_phasage_%s'%(tstr)+'.png')    



        

            




