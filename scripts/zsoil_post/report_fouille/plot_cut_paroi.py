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
import os,re


from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools
from zsoil_tools import postpro_lib as pl


def plot_cut_paroi(prob,coupes,etapes,ketapes,pvpath,figsize=(14,10),figpath='.'):
    tx = [etapes.tvect[kk] for kk in ketapes]
    etstr = [etapes.names[kk] for kk in ketapes]

    figpaths = []
    for k in range(len(coupes.names)):
        plane = vtk.vtkPlane()
        plane.SetOrigin(coupes.origins[k])
        plane.SetNormal(coupes.normals[k])
        loc_syst = np.array([np.cross(coupes.normals[k],(0,1,0)),(0,1,0),coupes.normals[k]])


        fig = plt.figure(figsize=figsize)
        bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
        fig.text(0.01,0.01,prob,size=8)

        AX = [fig.add_subplot(1,3,kk+1) for kk in range(3)]

        cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for kt in range(len(tx)):
            t = tx[kt]
            tstr = vtktools.get_tstr(t)
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(pvpath+'/'+prob+'_'+tstr+'_shell.vtu')
            reader.Update()
            grid = reader.GetOutput()

            segments = vtktools.get_section(grid,plane,matlist=[1,9,25,30],disp=True)
##                w = vtk.vtkPolyDataWriter()
##                w.SetFileName('cut_%s.vtk'%(titstr[k]))
##                w.SetInputData(grid)
##                w.Write()

            # disp:
            lab = False
            for seg in segments:
                if abs(seg[0][0])<1:
                    dy = 0.5*(seg[1][1]-seg[0][1])
                    if lab:
                        AX[0].plot([seg[2][3][0][2]*1e3,seg[2][3][1][2]*1e3],[seg[0][1],seg[1][1]],
                                   color=cols[kt])
                        AX[1].plot([seg[2][0][0]-dy*seg[2][2][0],seg[2][0][0]+dy*seg[2][2][0]],[seg[0][1],seg[1][1]],
                                   color=cols[kt])
                        AX[2].plot([0,seg[2][2][0],seg[2][2][0],0],[seg[0][1],seg[0][1],seg[1][1],seg[1][1]],
                                   color=cols[kt])
                    else:
                        AX[0].plot([seg[2][3][0][2]*1e3,seg[2][3][1][2]*1e3],[seg[0][1],seg[1][1]],
                                   color=cols[kt],label=etstr[kt])
                        AX[1].plot([seg[2][0][0]-dy*seg[2][2][0],seg[2][0][0]+dy*seg[2][2][0]],[seg[0][1],seg[1][1]],
                                   color=cols[kt],label=etstr[kt])
                        AX[2].plot([0,seg[2][2][0],seg[2][2][0],0],[seg[0][1],seg[0][1],seg[1][1],seg[1][1]],
                                   color=cols[kt],label=etstr[kt])
                        lab = True

        AX[0].set_xlabel(u'Déplacement horizontal [mm]')
        AX[1].set_xlabel(u'Moment [kNm/m]')
        AX[2].set_xlabel(u'Effort tranchant [kN/m]')
        for ax in AX:
            ax.set_ylabel(u'Altitude [msm]')
    ##                    ax.annotate('Max.: %1.0f mm'%(max(val1)),
    ##                                xy=(8,530),bbox=dict(fc='w',ec='k'))
    ##                    ax.annotate('Min.: %1.0f mm'%(min(val1)),
    ##                                xy=(65,530),bbox=dict(fc='w',ec='k'))
            ax.grid('on')
        AX[0].legend()
            
        fig.text(0.5,0.96,'Coupe '+coupes.names[k],size=16,ha='center')
            
        fig.tight_layout()
        fig.savefig(figpath+'/'+prob+'_dispParoi_'+re.sub(' ','_',coupes.names[k]))
        plt.close(fig)
        figpaths.append(figpath+'/'+prob+'_dispParoi_'+re.sub(' ','_',coupes.names[k]))

    return figpaths
    



        

            




