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
import matplotlib.tri as mtri
from matplotlib.ticker import FormatStrFormatter

import vtk
import os,re


from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools
from zsoil_tools import postpro_lib as pl

def plot_disp_cut_contours(prob,coupes,etapes,ketapes,pvpath,figsize=(14,6),figpath='.'):
    tx = [etapes.tvect[kk] for kk in ketapes]
    etstr = [etapes.names[kk] for kk in ketapes]

    figpaths = []
    
    for kt in range(len(tx)):
        t = tx[kt]
        tstr = vtktools.get_tstr(t)
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(pvpath+'/'+prob+'_'+tstr+'_vol.vtu')
        reader.Update()
        grid = reader.GetOutput()


        for k in range(len(coupes.names)):
            plane = vtk.vtkPlane()
            plane.SetOrigin(coupes.origins[k])
            plane.SetNormal(coupes.normals[k])
            loc_syst = np.array([np.cross(coupes.normals[k],(0,1,0)),(0,1,0)])

            val_y,crd,output = vtktools.get_section_vol(grid,plane,loc_syst,
                                                  array='DISP_TRA',component=1)
##            if coupes.normals[k][0]==0:
            val_x,crd,output = vtktools.get_section_vol(grid,plane,loc_syst,
                                                        array='DISP_TRA',component=0)
            val_z,crd,output = vtktools.get_section_vol(grid,plane,loc_syst,
                                                        array='DISP_TRA',component=2)
##            else:
##                val_x,crd,output = vtktools.get_section_vol(grid,plane,loc_syst,
##                                                      array='DISP_TRA',component=2)
##            xx = [pl.project_on_plane(loc_syst,coupes.origins[k],(crd[0][kk],crd[1][kk],crd[2][kk])) for kk in range(len(crd[0]
            triangles = [[output.GetCell(kc).GetPointIds().GetId(kk) for kk in range(3)] for kc in range(output.GetNumberOfCells())]
            triang = mtri.Triangulation(crd[0],crd[1],triangles)

##            w = vtk.vtkPolyDataWriter()
##            w.SetFileName('cut_%i.vtk'%(k))
##            w.SetInputData(output)
##            w.Write()

            fig = plt.figure(figsize=figsize)
            bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
            fig.text(0.01,0.01,prob,size=8)

            for kk in range(2):
                ax = fig.add_subplot(1,2,kk+1)
                if kk==0:
##                    val1 = [v*1e3 for v in val_x]
                    val1 = [vtktools.project_on_plane(loc_syst,[0,0,0],(val_x[kpt]*1e3,0,val_z[kpt]*1e3))[0] for kpt in range(len(val_y))]
                else:
                    val1 = [min(0,v)*1e3 for v in val_y]
##                    val1 = [v*1e3 for v in val_y]

                minmax = [1e10,-1e10]
                for kpt in range(len(crd[0])):
                    if crd[0][kpt]>coupes.bounds[k][0][0] and crd[0][kpt]<coupes.bounds[k][0][1]:
                        if crd[1][kpt]>coupes.bounds[k][1][0] and crd[1][kpt]<coupes.bounds[k][1][1]:
                            minmax[0] = min(minmax[0],val1[kpt])
                            minmax[1] = max(minmax[1],val1[kpt])

                cmap,norm,ticks = pl.GetDiscreteColormap(minmax)
                h = ax.tricontourf(triang,val1,cmap=cmap,norm=norm,levels=ticks,extend='both')

                cb = fig.colorbar(h)#,boundaries=ticks)
                if kk==0:
                    cb.set_label(u'Déplacement horizontal [mm]',size=14)
                    ax.annotate('Max.: %1.0f mm'%(max(val1)),
                                xy=(8,530),bbox=dict(fc='w',ec='k'))
                    ax.annotate('Min.: %1.0f mm'%(min(val1)),
                                xy=(65,530),bbox=dict(fc='w',ec='k'))
                else:
                    cb.set_label('Tassement [mm]',size=12)

                cb.ax.tick_params(labelsize=12)

                a0 = coupes.bounds[k][0][0]
                a1 = coupes.bounds[k][0][1]
                pt0 = [coupes.origins[k][0]+a0*loc_syst[0][0],coupes.origins[k][1]+a0*loc_syst[0][1],coupes.origins[k][2]+a0*loc_syst[0][2]]
                pt1 = [coupes.origins[k][0]+a1*loc_syst[0][0],coupes.origins[k][1]+a1*loc_syst[0][1],coupes.origins[k][2]+a1*loc_syst[0][2]]
                crd2D0 = vtktools.project_on_plane(loc_syst,coupes.origins[k],pt0)
                crd2D1 = vtktools.project_on_plane(loc_syst,coupes.origins[k],pt1)

                ax.set_xlim(crd2D0[0],crd2D1[0])
                ax.set_ylim(coupes.bounds[k][1][0],coupes.bounds[k][1][1])
                ax.set_ylabel(u'Altitude [msm]',size=12)
                ax.set_aspect('equal')
##                ax.annotate('Laupenstrasse',xy=(0.1,0.94),xycoords='axes fraction',size=12)
##                ax.annotate('Geleise',xy=(0.7,0.94),xycoords='axes fraction',size=12)
                ax.grid('both')

                
                fig.text(0.5,0.95,'Coupe '+coupes.names[k]+', %s'%(etstr[kt]),size=20,ha='center')
                
            fig.tight_layout()
            fig.savefig(figpath+'/'+prob+'_dispCut_'+re.sub(' ','_',coupes.names[k])+'_%1.0f'%(t))
            plt.close(fig)

            figpaths.append(figpath+'/'+prob+'_dispCut_'+re.sub(' ','_',coupes.names[k])+'_%1.0f'%(t))

    return figpaths
    



        

            




