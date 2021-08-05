# -*- coding: cp1252 -*-
# @description Plotting nodal results on cut section using vtk.
# @input zsoil results
# @output png
# @author Matthias Preisig
# @date 2018/06/22
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import vtk
import os


from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools
from zsoil_tools import postpro_lib as pl

pathname = '..'
pblist = ['M1224_3Dcomplet_v1_10']
tvect = [11]
tx = []

for kp,prob in enumerate(pblist):
    try:
        f = open(pathname+'/pv/'+prob+'.pvd')
        kt = 0
        for line in f:
##            if 'vol.vtu' in line:
            if '_vo_' in line:
                t = float(line[line.find('timestep=')+9:].split('"')[1])
                if t in tvect:
                    tx.append(t)
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
                res.read_s01()  # volumics
        ptemp = os.getcwd()
        os.chdir('../pv')
        vtktools.write_vtu(res,vol=True,verbose=False)
        os.chdir(ptemp)

    figname = prob
    
    for kt in range(len(tx)):
        t = tx[kt]
##        tstr = str(int(t)).rjust(3,'0')+'_'+str(int((t-int(t))*100)).rjust(2,'0')
        tstr = str(int(t))+'_'+str(int((t-int(t))*100))
        reader = vtk.vtkXMLUnstructuredGridReader()
##        reader.SetFileName(pathname+'/pv/'+prob+'_'+tstr+'_vol.vtu')
        reader.SetFileName(pathname+'/pv/'+prob+'_vo_'+tstr+'.vtu')
        reader.Update()
        grid = reader.GetOutput()


        ##############################################################################
        # coupes:
        ##############################################################################
        titstr = ['3-3','4-4','4\'-4\'','B11-B12']
        orig = [[26698.86920,0,-11665.01342],
                [26716.60017,0,-11677.80772],
                [26707.39172,0,-11671.22828],
                [2.67499e+04,0,-1.16753e+04]]
        normal = [[0.78297519,0,-0.62205293],
                  [0.82629588,0,-0.56323628],
                  [0.82629588,0,-0.56323628],
                  [0.89442719,0,0.4472136]]
        BOUNDS = [[-50,40],
                  [-60,40],
                  [-55,45],
                  [-60,40]]

##        fig = plt.figure(figsize=(14,16))
##        bbox_props = dict(boxstyle="round", fc="w", ec="k", lw=2)
##        fig.text(0.01,0.01,prob,size=8)
##
##        ax = plt.subplot(111)
##        for ke in range(res.nShells):
##            inel = res.shell.inel[ke]
##            ax.plot([res.coords[0][kn-1] for kn in inel],[-res.coords[2][kn-1] for kn in inel],'k')
##
##        for k in range(len(orig)):
##            a0 = BOUNDS[k][0]
##            a1 = BOUNDS[k][1]
##            loc_syst = numpy.array([numpy.cross(normal[k],(0,1,0)),(0,1,0)])
##            ax.plot([orig[k][0]+a0*loc_syst[0][0],orig[k][0]+a1*loc_syst[0][0]],
##                    [-orig[k][2]-a0*loc_syst[0][2],-orig[k][2]-a1*loc_syst[0][2]],'r')
##        
##        fig.tight_layout()
##        fig.savefig(prob+'_layout')
##        plt.close(fig)

        for k in range(len(orig)):
            plane = vtk.vtkPlane()
            plane.SetOrigin(orig[k])
            plane.SetNormal(normal[k])
            loc_syst = np.array([np.cross(normal[k],(0,1,0)),(0,1,0)])

            val_y,crd,output = vtktools.get_section_vol(grid,plane,loc_syst,
                                                  array='disp',component=1)
            if normal[k][0]==0:
                val_x,crd,output = vtktools.get_section_vol(grid,plane,loc_syst,
                                                      array='disp',component=0)
            else:
                val_x,crd,output = vtktools.get_section_vol(grid,plane,loc_syst,
                                                      array='disp',component=2)
            w = vtk.vtkPolyDataWriter()
            w.SetFileName('cut_%i.vtk'%(k))
            w.SetInputData(output)
            w.Write()

            fig = plt.figure(figsize=(9,8))
            bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
            fig.text(0.01,0.01,prob,size=8)

            for kk in range(2):
                ax = fig.add_subplot(2,1,kk+1)
                if kk==0:
                    val1 = [-v*1e3 for v in val_x]
                else:
                    val1 = [min(0,v)*1e3 for v in val_y]

                CS = pl.contourf(ax,val1,crd,output,loc_syst,orig[k])
                
                cb = fig.colorbar(CS)
                cbar_min = cb.get_clim()[0]
                cbar_max = cb.get_clim()[1]
                cbar_step = (cbar_max-cbar_min)/(len(CS.levels)-1)
                cb.ax.set_yticklabels(['{:.1f}'.format(x) for x in np.arange(cbar_min, cbar_max+cbar_step, cbar_step)])#, fontsize=16, weight='bold')
                if kk==0:
                    cb.set_label(u'Déplacement horizontal [mm]',size=14)
                    ax.annotate('Max.: %1.0f mm'%(max(val1)),
                                xy=(8,530),bbox=dict(fc='w',ec='k'))
                    ax.annotate('Min.: %1.0f mm'%(min(val1)),
                                xy=(65,530),bbox=dict(fc='w',ec='k'))
                else:
                    cb.set_label('Tassement [mm]',size=12)
##                    ax.annotate('Min.: %1.0f mm'%(min([val1[kk] for kk in range(len(crd[0])) if crd[0][kk]>40])),
##                                xy=(65,530),bbox=dict(fc='w',ec='k'))
##                    ax.annotate('Min.: %1.0f mm'%(min([val1[kk] for kk in range(len(crd[0])) if crd[0][kk]<40])),
##                                xy=(8,530),bbox=dict(fc='w',ec='k'))

                cb.ax.tick_params(labelsize=12)

                a0 = BOUNDS[k][0]
                a1 = BOUNDS[k][1]
                pt0 = [orig[k][0]+a0*loc_syst[0][0],orig[k][1]+a0*loc_syst[0][1],orig[k][2]+a0*loc_syst[0][2]]
                pt1 = [orig[k][0]+a1*loc_syst[0][0],orig[k][1]+a1*loc_syst[0][1],orig[k][2]+a1*loc_syst[0][2]]
                crd2D0 = vtktools.project_on_plane(loc_syst,orig[k],pt0)
                crd2D1 = vtktools.project_on_plane(loc_syst,orig[k],pt1)

                ax.set_xlim(crd2D0[0],crd2D1[0])
                ax.set_ylim(min(crd[1]),max(crd[1]))
                ax.set_ylabel(u'Altitude [msm]',size=12)
                ax.set_aspect('equal')
##                ax.annotate('Laupenstrasse',xy=(0.1,0.94),xycoords='axes fraction',size=12)
##                ax.annotate('Geleise',xy=(0.7,0.94),xycoords='axes fraction',size=12)

                
                ax.set_title('Coupe '+titstr[k],size=20)
                
            fig.tight_layout()
            fig.savefig(prob+'_dispCut_'+titstr[k].split()[0])
            plt.close(fig)

    



        

            




