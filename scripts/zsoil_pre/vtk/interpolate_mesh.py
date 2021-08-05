import vtk
import numpy as np
import math

layers = ['Couche 7d',
          'Couche 4a-4b',
          'Couche 3c',
          'TN HKD']

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

def getNewYCoordOfLayer(x,y0,z,kl):      
    rayStart = [x,y0-100,z]
    rayEnd = [x,y0+100,z]
    t = vtk.mutable(0)
    xyz = [0.0, 0.0, 0.0]
    pcoords = [0.0, 0.0, 0.0]
    subId = vtk.mutable(0)
    if locators[kl].IntersectWithLine(rayStart,rayEnd,0.0001,t,xyz,pcoords,subId):
        return xyz[1]
    else:
##        print 'point (%1.1f;%1.1f) not inside domain'%(x,z)
        return y0


f = open('../recu/troinexTest_boreholes.inp')
of = open('troinexTest_boreholes_1.inp','w')

ynew = range(10,100)
# vertical nodal spacing in original mesh:
dy = 0.5

of.write(f.next())
line = f.next()
of.write(line)
nNodes = int(line.split()[5])
for line in f:
    if '.ing' in line:
        of.write(line)
        for kn in range(nNodes):
            line = f.next()
            v = line.split()
            x = float(v[1])
            y0 = float(v[2])
            z = float(v[3])
            # index of node layer:
            yindex = int((y0-414)/dy)
            if yindex<=4:   # 4 elements in lowest layer
                ybot = getNewYCoordOfLayer(x,y0,z,0)
                ytop = getNewYCoordOfLayer(x,y0,z,1)
                y1 = ybot + yindex/4.*(ytop-ybot)
            elif yindex<=6: # 2 elements in 2nd layer
                ybot = getNewYCoordOfLayer(x,y0,z,1)
                ytop = getNewYCoordOfLayer(x,y0,z,2)
                y1 = ybot + (yindex-4)/2.*(ytop-ybot)
            else:   # 8 elements in top layer
                ybot = getNewYCoordOfLayer(x,y0,z,2)
                ytop = getNewYCoordOfLayer(x,y0,z,3)
                y1 = ybot + (yindex-6)/8.*(ytop-ybot)
            of.write('%i %s %1.12f %s 0\n'%(kn+1,v[1],y1,v[3]))
            if kn%10000==9999:
                print '%i nodes out of %i'%(kn,nNodes)
    else:
        of.write(line)
f.close()
of.close()
