import vtk
import numpy as np
import math
from zsoil_tools import zsoil_inp as zi

pathname = '.'
prob = 'M1224_3Dcomplet_v1_6'
mesh = zi(pathname,prob)
mesh.read_inp()

layers = ['remb_CHUV',
          'colluvions',
          'molasse']

locators = []
for kl,layer in enumerate(layers):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(layer+'.vtp')
    reader.Update()
    output = reader.GetOutput()
    cellLocator = vtk.vtkCellLocator()
    cellLocator.SetDataSet(output)
    cellLocator.BuildLocator()
    locators.append(cellLocator)

def getLayerRelPos(x,y,z,kl):      
    rayStart = [x,y-100,z]
    rayEnd = [x,y+100,z]
    t = vtk.mutable(0)
    xyz = [0.0, 0.0, 0.0]
    pcoords = [0.0, 0.0, 0.0]
    subId = vtk.mutable(0)
    if locators[kl].IntersectWithLine(rayStart,rayEnd,0.0001,t,xyz,pcoords,subId):
        if xyz[1]>y:
            return -1
        elif xyz[1]<=y:
            return 1
        else:
            return 0
    else:
##        print 'point (%1.1f;%1.1f) not inside domain'%(x,z)
        return 0


f = open('M1224_3Dcomplet_v1_6.inp')
of = open('M1224_3Dcomplet_v1_7.inp','w')

of.write(next(f))
line = next(f)
of.write(line)
for line in f:
    if '.i0g' in line:
        of.write(line)
        for ke in range(mesh.nVolumics):
            xm = np.mean([mesh.coords[0][kn-1] for kn in mesh.vol.inel[ke]])
            ym = np.mean([mesh.coords[1][kn-1] for kn in mesh.vol.inel[ke]])
            zm = np.mean([mesh.coords[2][kn-1] for kn in mesh.vol.inel[ke]])
            l0 = getLayerRelPos(xm,ym,zm,0) # remb CHUV
            l1 = getLayerRelPos(xm,ym,zm,1) # colluvions
            l2 = getLayerRelPos(xm,ym,zm,2) # molasse
            if l2==-1:
                mat = 5
            elif l1==1:
                if l0==-1:
                    mat = 3
                else:
                    mat = 2
            elif l1==-1 and l2==1:
                mat = 4
            else:
                mat = 5

            line = next(f)
            v = line.split()
            astr = ''
            for kk in range(14):
                astr += v[kk]+' '
            astr += str(mat)
            for kk in range(15,20):
                astr += ' '+v[kk]
            astr += '\n'
            of.write(astr)
    else:
        of.write(line)
f.close()
of.close()
