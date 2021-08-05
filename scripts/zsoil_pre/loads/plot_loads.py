# -*- coding: cp1252 -*-
import sys,math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

from zsoil_tools import zsoil_inp as zi


pathname = '//192.168.1.51/Mandats sur H RAID0/M1013_Hardturm/v1017/West'
prob = 'M1013_W_v6_pilesV1017_cent1_l161213=844MN'
pathname = '//192.168.1.51/Mandats sur H RAID0/M1013_Hardturm/v1017/Ost'
prob = 'M1013_O_v8_pilesV1017_25m90_l170220=871MN'

mesh = zi(pathname,prob)
mesh.read_inp()

bbox_props = dict(boxstyle="round,pad=0", fc="w", ec="k", lw=0.4)
fig = plt.figure(figsize=(8,10))
fig.text(0.01,0.01,prob,size=8)

ax = fig.add_subplot(111)

time = 5

patches = []
cvect = []
scale = 1e-4
scale2 = 5e-6
nodal = 0
nLx = 0
nLz = 0
for kn in range(mesh.nNodalLoads):
    load = mesh.nodalLoads[kn][1][1]*mesh.LFs[mesh.nodalLoads[kn][2]].at(time)
    nodal += load
    x = mesh.coords[0][mesh.nodalLoads[kn][0]-1]
    z = mesh.coords[2][mesh.nodalLoads[kn][0]-1]
    nLx += load*x
    nLz += load*z
    patches.append(Circle((x,-z),radius=load*scale,fc='b'))
    cvect.append(0)
cgNodal = (nLx/nodal,-nLz/nodal)
patches.append(Circle(cgNodal,radius=nodal*scale2,fc='b'))
cvect.append(1)

beam = 0
bLx = 0
bLz = 0
scbl = 2e-4
for kb in range(mesh.nBeamLoads):
    ke = mesh.beamLoads[kb][0]
    inel = mesh.beam.inel[ke-1]
    x0 = mesh.coords[0][inel[0]-1]
    x1 = mesh.coords[0][inel[1]-1]
    z0 = -mesh.coords[2][inel[0]-1]
    z1 = -mesh.coords[2][inel[1]-1]
    xm = 0.5*(x0+x1)
    zm = 0.5*(z0+z1)
    dx = x0-x1
    dz = z0-z1
    dl = (dx**2+dz**2)**0.5
    l0 = mesh.beamLoads[kb][1][1]*mesh.LFs[mesh.beamLoads[kn][2]].at(time)
    l1 = mesh.beamLoads[kb][1][4]*mesh.LFs[mesh.beamLoads[kn][2]].at(time)
    load = dl*0.5*(l0+l1)
    beam += load
    bLx += load*xm
    bLz += load*zm
    if abs(dx)<1e-2:
        ax.plot([x0,x1,x1+l1*scbl,x0+l0*scbl,x0],
                [z0,z1,z1,z0,z0],'k',linewidth=0.5)
    elif abs(dz)<1e-2:
        ax.plot([x0,x1,x1,x0,x0],
                [z0,z1,z1+l1*scbl,z0+l0*scbl,z0],'k',linewidth=0.5)
    else:
        ax.plot([x0,x1],
                [z0,z1],'k',linewidth=0.5)
        print('Lasten nicht orthogonal! (%1.2e,%1.2e)'%(dx,dz))
cgBeam = (bLx/beam,bLz/beam)

patches.append(Circle(cgBeam,radius=beam*scale2,fc='r'))
cvect.append(2)
cmap = plt.cm.get_cmap('jet')
pc = PatchCollection(patches,edgecolors=('k',),alpha=0.8,cmap=cmap)
pc.set_array(np.array(cvect))
ax.add_collection(pc)

ax.annotate(u'Summe Stützenlasten: %1.1f MN'%(nodal*1e-3),
            xy=cgNodal,ha='center',size=14)
ax.annotate('Summe Linienlasten: %1.1f MN'%(beam*1e-3),
            xy=cgBeam,ha='center',size=14)

ax.set_aspect('equal')

fig.tight_layout()
fig.savefig(prob+'_lastverteilung')
plt.close(fig)
