# -*- coding: cp1252 -*-
import sys,math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from zsoil_tools import zsoil_inp as zi


pathname = '..'
prob = 'M1279_Scheibe_DN600_v3'
##pathname = '../../scheibe_v2'
##prob = 'M1279_Scheibe_DN600_v2_2'

mesh = zi(pathname,prob)
mesh.read_inp(sections=['ing','ibg','ilg','ics','icg','inb',
                        'pob','gob','inl','ibf','gsl','itg',
                        'brc','ipg','i0g'])

plot_beamLoads = False
plot_surfLoads = True
plot_nodalLoads = True

bbox_props = dict(boxstyle="round,pad=0", fc="w", ec="k", lw=0.4)
fig = plt.figure(figsize=(4,20))
fig.text(0.01,0.01,prob,size=8)

ax = fig.add_subplot(111)

time = 5

patches = []
cvect = []
scale = 1e-4
scale2 = 5e-6

if plot_nodalLoads:
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
        ax.annotate('%1.0f kN'%(mesh.nodalLoads[kn][1][1]),
                    xy=(x,-z),size=8)
    ##    ax.annotate('%s: %1.0f kN'%(mesh.nodalLoads[kn][3].split()[1],mesh.nodalLoads[kn][1][1]),
    ##                xy=(x,-z))
        cvect.append(0)
    cgNodal = (nLx/nodal,-nLz/nodal)
    patches.append(Circle(cgNodal,radius=nodal*scale2,fc='b'))
    cvect.append(1)
    patches.append(Circle((188,290),radius=nodal*scale2,fc='b'))
    cvect.append(1)

if plot_beamLoads:
    beam = 0
    bLx = 0
    bLz = 0
    scbl = 1e-3
    bloads = []
    bload_pos = []
    bload_labels = []
    bload_dir = []
    for kb in range(mesh.nBeamLoads):
        ke = mesh.beamLoads[kb][0]
        inel = mesh.beam.inel[ke-1]
        y0 = mesh.coords[1][inel[0]-1]
        x0 = mesh.coords[0][inel[0]-1]
        x1 = mesh.coords[0][inel[1]-1]
        z0 = -mesh.coords[2][inel[0]-1]
        z1 = -mesh.coords[2][inel[1]-1]
        xm = 0.5*(x0+x1)
        zm = 0.5*(z0+z1)
        dx = x0-x1
        dz = z0-z1
        dl = (dx**2+dz**2)**0.5

        l0 = mesh.beamLoads[kb][1][1]*mesh.LFs[mesh.beamLoads[kb][2]].at(time)
        l1 = mesh.beamLoads[kb][1][4]*mesh.LFs[mesh.beamLoads[kb][2]].at(time)
        load = dl*0.5*(l0+l1)
        beam += load
        bLx += load*xm
        bLz += load*zm
        n = np.array([-dz,dx])/dl
        ax.plot([x0,x1,x1+n[0]*l1*scbl,x0+n[0]*l0*scbl,x0],
                [z0,z1,z1+n[1]*l1*scbl,z0+n[1]*l0*scbl,z0],'k',linewidth=0.5)
        if not mesh.beamLoads[kb][3] in bload_labels:
            bload_labels.append(mesh.beamLoads[kb][3])
            bload_pos.append([np.array([xm,zm])])
            bloads.append(load)
            bload_dir.append(x1-x0)
        else:
            ind = bload_labels.index(mesh.beamLoads[kb][3])
            bload_pos[ind].append(np.array([xm,zm]))
            bloads[ind] += load
    cgBeam = np.array([bLx/beam,bLz/beam])

    for kb,bload in enumerate(bloads):
        if abs(bload_dir[kb])<1e-3:
            ax.annotate('%1.0f kN'%(bload),
                        xy=np.mean(bload_pos[kb],0),va='center',rotation=90)
        else:
            ax.annotate('%1.0f kN'%(bload),
                        xy=np.mean(bload_pos[kb],0),ha='center',rotation=0)

    patches.append(Circle(cgBeam,radius=beam*scale2,fc='r'))
    cvect.append(2)
    patches.append(Circle((3,-32),radius=beam*scale2,fc='r'))
    cvect.append(2)

if plot_surfLoads:
    surf = 0
    surf0 = 0
    sLx = 0
    sLz = 0
    for ks in range(mesh.nSurfaceLoads):
        mesh.surfLoads[ks].compute_resultant(mesh)
        load = mesh.surfLoads[ks].resultant[1]*mesh.LFs[mesh.surfLoads[ks].LF].at(time)

        if mesh.surfLoads[ks].LF==1:
            surf += load
            sLx += load*mesh.surfLoads[ks].CG[0]
            sLz += load*mesh.surfLoads[ks].CG[2]
            for f in mesh.surfLoads[ks].faces:
                patches.append(Polygon([(mesh.coords[0][f[1][k]-1],-mesh.coords[2][f[1][k]-1]) for k in range(4)]))
                cvect.append(4)
            ax.annotate('%1.0f kN'%(load),
                        xy=(mesh.surfLoads[ks].CG[0],-mesh.surfLoads[ks].CG[2]),bbox=bbox_props,
                        ha='center')
        else:
            surf0 += load
    print(surf0)
    cgSurf = np.array([sLx/surf,-sLz/surf])

cmap = plt.cm.get_cmap('jet')
pc = PatchCollection(patches,edgecolors=('k',),alpha=0.8,cmap=cmap)
pc.set_array(np.array(cvect))
ax.add_collection(pc)

ax.annotate(u'Schwerpunkt Stützenlasten ($\sum=$%1.1f MN)'%(nodal*1e-3),
            xy=(3,-110),ha='left',size=12,rotation=90)
ax.annotate(u'Schwerpunkt Deckenlasten ($\sum=$%1.1f MN)'%(surf*1e-3),
            xy=(4.5,-110),ha='left',size=12,rotation=90)
##ax.annotate('Schwerpunkt Wand+Kernlasten ($\sum=$%1.1f MN)'%((beam)*1e-3),
##            xy=(190,287),ha='left',size=12)
ax.annotate(u'Total Lasten: %1.1f MN'%((nodal+surf)*1e-3),
            xy=(6,-110),ha='left',size=12,rotation=90)

ax.set_aspect('equal')
ax.set_xlim([-6,6])
ax.set_ylim([-127,1])
ax.grid('on')

fig.tight_layout()
fig.savefig(prob+'_lastverteilung')
plt.close(fig)
