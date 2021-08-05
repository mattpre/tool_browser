# -*- coding: cp1252 -*-
import numpy,os,math,cmath
from numpy import linalg as la
from numpy import fft
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import cPickle as pickle
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from openpyxl import load_workbook
from scipy import interpolate
import six
from matplotlib import colors
import vtk

import sys
##sys.path.append('D:/Mandats/python/zsoil_tools')
##sys.path.append('//JAN-09/Mandats/Mandats/python/zsoil_tools')
from zsoil_tools import zsoil_results as zr
import time


pathname = r'\\192.168.1.55\Mandats\M1086 Tridel\calc'
pblist = ['tridel_v7_demi']
pblist = ['tridel_v7_demi_c100kPaMA']
pblist = ['tridel_v5_demi_cnt_shellhinge_Levas_Radier_domaine_C3Dbeam']


for kf,prob in enumerate(pblist):
    figname = prob
    try:
        res = pickle.load(open(prob+'_evol.p', "rb" ))
        print prob+' loaded'
        tsteps = res.out_steps
        tvect = [res.steps[kt].time for kt in tsteps]
    except:
        res = zr(pathname,prob)
        res.read_rcf()
        res.read_his()

        tvect = [3,4,5,6,7,8,9,10,11,12]
        tsteps = []
        for kt,step in enumerate(res.steps):
            if step.time in tvect and step.sf==0 and step.conv_status==-1:
                tsteps.append(kt)
##        tsteps = [5,6]
                
        res.out_steps = tsteps
        res.read_dat()
        res.read_s04()
        res.read_s00()
        pickle.dump(res, open(prob+'_evol.p', "wb" ) ,pickle.HIGHEST_PROTOCOL)

    nails = []
    Nmax = []
    Nmin = []
    for kn in range(res.nNails):
        nail = res.nails[kn]
        Nmax.append([])
        Nmin.append([])
        nails.append([])
        for kb in nail.beams:
            nails[-1].append(res.num_beams.index(kb))
            Nmax[-1].append(-1e10)
            Nmin[-1].append(1e10)

    for kkt,kt in enumerate(tsteps):
        step = res.steps[kt]
        for kn in range(res.nNails):
            for kkb,kb in enumerate(nails[kn]):
                Nmax[kn][kkb] = max(Nmax[kn][kkb],step.beam.force[0][kb])
                Nmin[kn][kkb] = min(Nmin[kn][kkb],step.beam.force[0][kb])
    
if False:
    pd = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    arr = vtk.vtkFloatArray()
    arr.SetName('Nmax')

    scale = 0.005
    kpt = 0
    for kn,nail in enumerate(nails):
        for kkb,kb in enumerate(nail):
            inel = res.beam.inel[kb]
            arr.InsertNextTuple1(Nmax[kn][kkb])
            cell = vtk.vtkQuad()
            Ids = cell.GetPointIds()
            Ids.SetId(0,kpt)
            Ids.SetId(1,kpt+1)
            Ids.SetId(2,kpt+2)
            Ids.SetId(3,kpt+3)
            kpt += 4
            points.InsertNextPoint(res.coords[0][inel[0]-1],
                                   res.coords[1][inel[0]-1],
                                   res.coords[2][inel[0]-1])
            points.InsertNextPoint(res.coords[0][inel[0]-1]+Nmax[kn][kkb]*scale,
                                   res.coords[1][inel[0]-1],
                                   res.coords[2][inel[0]-1])
            points.InsertNextPoint(res.coords[0][inel[1]-1]+Nmax[kn][kkb]*scale,
                                   res.coords[1][inel[1]-1],
                                   res.coords[2][inel[1]-1])
            points.InsertNextPoint(res.coords[0][inel[1]-1],
                                   res.coords[1][inel[1]-1],
                                   res.coords[2][inel[1]-1])
            cells.InsertNextCell(cell)

    pd.SetPoints(points)
    pd.SetPolys(cells)
    pd.GetCellData().AddArray(arr)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName('Nmax.vtp')
    writer.SetInputData(pd)

    writer.Write()
    pd = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    arr = vtk.vtkFloatArray()
    arr.SetName('Nmax')

    scale = 0.005
    kpt = 0
    for kn,nail in enumerate(nails):
        for kkb,kb in enumerate(nail):
            inel = res.beam.inel[kb]
            arr.InsertNextTuple1(Nmin[kn][kkb])
            cell = vtk.vtkQuad()
            Ids = cell.GetPointIds()
            Ids.SetId(0,kpt)
            Ids.SetId(1,kpt+1)
            Ids.SetId(2,kpt+2)
            Ids.SetId(3,kpt+3)
            kpt += 4
            points.InsertNextPoint(res.coords[0][inel[0]-1],
                                   res.coords[1][inel[0]-1],
                                   res.coords[2][inel[0]-1])
            points.InsertNextPoint(res.coords[0][inel[0]-1]+Nmin[kn][kkb]*scale,
                                   res.coords[1][inel[0]-1],
                                   res.coords[2][inel[0]-1])
            points.InsertNextPoint(res.coords[0][inel[1]-1]+Nmin[kn][kkb]*scale,
                                   res.coords[1][inel[1]-1],
                                   res.coords[2][inel[1]-1])
            points.InsertNextPoint(res.coords[0][inel[1]-1],
                                   res.coords[1][inel[1]-1],
                                   res.coords[2][inel[1]-1])
            cells.InsertNextCell(cell)

    pd.SetPoints(points)
    pd.SetPolys(cells)
    pd.GetCellData().AddArray(arr)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName('Nmin.vtp')
    writer.SetInputData(pd)

    writer.Write()

if True:
    fig = plt.figure(figsize=(16,12))
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.4)
    fig.text(0.01,0.01,prob,size=8)

    ax = fig.add_subplot(111)

    lstr = {5:'b.-',6:'r.-'}
##    lstr = {13:'b.-',14:'r.-',15:'g.-'}
    scale = 0.002
    kpt = 0
    for kn,nail in enumerate(nails):
        y0 = res.coords[1][res.beam.inel[nail[0]][0]-1]
        z0 = res.coords[2][res.beam.inel[nail[0]][0]-1]
        if y0>5:
            y00 = 2*9
        elif y0>4:
            y00 = 9
        else:
            y00 = 0
        for kkb,kb in enumerate(nail):
            inel = res.beam.inel[kb]
            x = [res.coords[0][inel[0]-1]+Nmin[kn][kkb]*scale,
                 res.coords[0][inel[0]-1]+Nmax[kn][kkb]*scale,
                 res.coords[0][inel[1]-1]+Nmax[kn][kkb]*scale,
                 res.coords[0][inel[1]-1]+Nmin[kn][kkb]*scale,
                 res.coords[0][inel[0]-1]+Nmin[kn][kkb]*scale]
            y = [y00+((res.coords[1][inel[0]-1]-y0)**2+(res.coords[2][inel[0]-1]-z0)**2)**0.5,
                 y00+((res.coords[1][inel[0]-1]-y0)**2+(res.coords[2][inel[0]-1]-z0)**2)**0.5,
                 y00+((res.coords[1][inel[1]-1]-y0)**2+(res.coords[2][inel[1]-1]-z0)**2)**0.5,
                 y00+((res.coords[1][inel[1]-1]-y0)**2+(res.coords[2][inel[1]-1]-z0)**2)**0.5,
                 y00+((res.coords[1][inel[0]-1]-y0)**2+(res.coords[2][inel[0]-1]-z0)**2)**0.5]
            ax.plot(x,y,'k')
            ax.plot([res.coords[0][inel[0]-1],
                     res.coords[0][inel[0]-1]],[y00,y00+8],'k')
            ax.annotate('min: %1.0f kN\nmax: %1.0f kN'%(min(Nmin[kn]),max(Nmax[kn])),
                        xy=(res.coords[0][inel[0]-1],y00-1),
                        ha='center',size=12)


##    ax.legend(loc='lower right')
    ax.grid(b=True,which='both',axis='y')
##    ax.set_xlim([5,6])
##    ax.set_xlabel('Etappe')
##    ax.set_ylabel('Setzung Pfahlkopf [mm]')
    fig.tight_layout()
    fig.savefig(prob+'_clous_envel')
    plt.close(fig)

