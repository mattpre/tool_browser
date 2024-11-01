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
import os,re

from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

class coupe:
    def __init__(self,name,orig,norm):
        self.name = name
        self.origin = orig
        self.normal = norm
coupes = []
coupes.append(coupe('Côté mur aile 1',[5.25,0,4],[1,0,0]))
coupes.append(coupe('Côté mur aile 2',[6.5,0,4],[1,0,0]))


pathname = '../..'
pvpath = '../../pv'
pblist = ['M1396_Est_3D_v3_2c_ortho_excMurVoileC_ELU']
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
        res.read_dat()
        tsteps = []
        t0step = -1
        for kt,step in enumerate(res.steps):
            if step.conv_status==-1:
##                if step.time in etapes.tvect or abs(step.time%1)<1e-3:
                tsteps.append(kt)
                if step.time==4:
                    t0step = kt
        res.out_steps = tsteps
        res.read_s00()
        for lab in res.ele_group_labels:
            if lab=='SHELLS':
                res.read_s02()
        ptemp = os.getcwd()
        os.chdir(pvpath)
        pvdFile = vtktools.create_pvd(res)
        vtktools.write_vtu(res,shells=True,verbose=False,pvdFile=pvdFile,refstep=res.steps[t0step])
##        vtktools.write_vtu(res,vol=True,verbose=False,pvdFile=pvdFile)
        vtktools.save_pvd(pvdFile)
        os.chdir(ptemp)

    t = 17
    tstr = vtktools.get_tstr(t,t0=4)
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(pvpath+'/'+prob+'_'+tstr+'_shell.vtu')
    reader.Update()
    mesh = reader.GetOutput()

    for kc in range(len(coupes)):
        fig = plt.figure(figsize=(10,8))
        bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
        fig.text(0.01,0.01,prob,size=8)

        AX = [fig.add_subplot(1,3,kk+1) for kk in range(3)]
        plane = vtk.vtkPlane()
        plane.SetOrigin(coupes[kc].origin)
        plane.SetNormal(coupes[kc].normal)
            
        segments = vtktools.get_section(mesh,plane,matlist=[13],disp=True)
        lab = False
        kt = 0
        etstr=['aaa']
        for seg in segments:
            if abs(seg[0][0])<1:
                dy = 0.5*(seg[1][1]-seg[0][1])
                if lab:
                    AX[0].plot([seg[2][3][0][0]*1e3,seg[2][3][1][0]*1e3],[seg[0][1],seg[1][1]],
                               color=colors[kt])
                    AX[1].plot([seg[2][0][0]-dy*seg[2][2][0],seg[2][0][0]+dy*seg[2][2][0]],[seg[0][1],seg[1][1]],
                               color=colors[kt])
                    AX[2].plot([0,seg[2][2][0],seg[2][2][0],0],[seg[0][1],seg[0][1],seg[1][1],seg[1][1]],
                               color=colors[kt])
                else:
                    AX[0].plot([seg[2][3][0][0]*1e3,seg[2][3][1][0]*1e3],[seg[0][1],seg[1][1]],
                               color=colors[kt],label=etstr[kt])
                    AX[1].plot([seg[2][0][0]-dy*seg[2][2][0],seg[2][0][0]+dy*seg[2][2][0]],[seg[0][1],seg[1][1]],
                               color=colors[kt],label=etstr[kt])
                    AX[2].plot([0,seg[2][2][0],seg[2][2][0],0],[seg[0][1],seg[0][1],seg[1][1],seg[1][1]],
                               color=colors[kt],label=etstr[kt])
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
            ax.set_title(' ')
        AX[0].legend()
            
        fig.text(0.5,0.97,'Coupe '+coupes[kc].name,size=14,ha='center')
            
        fig.tight_layout()
        fig.savefig(prob+'_dispParoi_'+re.sub(' ','_',coupes[kc].name))
        plt.close(fig)

    



        

            




