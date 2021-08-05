import vtk
import numpy as np
import math

layers = ['terrain']

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

x0 = 0
z0 = 0

def getNewYCoordOfLayer(x,y0,z,kl=0):
    rayStart = [x+x0,y0-1000,z]
    rayEnd = [x+x0,y0+1000,z]
##    rayStart = [x,y0-1000,z]
##    rayEnd = [x,y0+1000,z]
    t = vtk.mutable(0)
    xyz = [0.0, 0.0, 0.0]
    pcoords = [0.0, 0.0, 0.0]
    subId = vtk.mutable(0)
##    print(x,y0,z)
##    print(locators[kl].IntersectWithLine(rayStart,rayEnd,0.0001,t,xyz,pcoords,subId))
    if locators[kl].IntersectWithLine(rayStart,rayEnd,0.0001,t,xyz,pcoords,subId):
##        print(rayStart)
        return xyz[1]
    else:
##        print 'point (%1.1f;%1.1f) not inside domain'%(x,z)
        return y0


f = open('M1224_3Dcomplet_v1.inp')
of = open('M1224_3Dcomplet_v1_1.inp','w')

# vertical nodal spacing in original mesh:

of.write(next(f))
line = next(f)
of.write(line)
nNodes = int(line.split()[5])
for line in f:
    if '.ing' in line:
        of.write(line)
        for kn in range(nNodes):
            line = next(f)
##            print(line)
            v = line.split()
            x = float(v[1])
            y0 = float(v[2])
            z = float(v[3])
            # index of node layer:
            if abs(y0-720)<0.01:
                y1 = getNewYCoordOfLayer(x,y0,z,0)
            else:
                y1 = y0
            of.write('%i %s %1.12f %s 0\n'%(kn+1,v[1],y1,v[3]))
            if kn%10000==9999:
                print('%i nodes out of %i'%(kn,nNodes))
    else:
        of.write(line)
f.close()
of.close()
