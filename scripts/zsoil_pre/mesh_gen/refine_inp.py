import vtk
import numpy as np
from zsoil_tools import zsoil_inp as zi

pathname = '..'
prob = 'M1224_3Dcomplet_v1_4'
femesh = zi(pathname,prob)
femesh.read_inp()


points = vtk.vtkPoints()
for kpt in range(femesh.nNodes):
    points.InsertNextPoint(femesh.coords[0][kpt],
                           femesh.coords[1][kpt],
                           femesh.coords[2][kpt])
ca0 = vtk.vtkCellArray()
for ke in range(femesh.nVolumics):
    hexa = vtk.vtkHexahedron()
    for kn in range(8):
        hexa.GetPointIds().SetId(kn,femesh.vol.inel[ke][kn]-1)
    ca0.InsertNextCell(hexa)
mesh0 = vtk.vtkUnstructuredGrid()
mesh0.SetPoints(points)
mesh0.SetCells(vtk.VTK_HEXAHEDRON,ca0)

pd = vtk.vtkPolyData()
pd.SetPoints(points)
pointLocator = vtk.vtkPointLocator()
pointLocator.SetDataSet(pd)
pointLocator.InitPointInsertion(points,points.GetBounds())
pointLocator.BuildLocator()

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('mesh0.vtu')
writer.SetInputData(mesh0)
writer.Write()

##facenodes = []
##f = open('face_nodes1.txt')
##for line in f:
##    v = line.split()
##    facenodes.append(int(float(v[0])))
##f.close()
##elist = []
##f = open('elements.txt')
##for line in f:
##    elist.append(int(float(line)))
##f.close()
for fg in femesh.face_groups:
    if fg.name=='master':
        break

iedict = {'0123v':[6,2,1,5,7,3,0,4],
          '0123':[7,3,2,6,4,0,1,5],
          '4567v':[1,5,6,2,0,4,7,3],
          '4567':[2,6,7,3,1,5,4,0],
          '0145v':[6,5,4,7,2,1,0,3],
          '0145':[7,4,0,3,6,5,1,2],
          '2367v':[1,2,3,0,5,6,7,4],
          '2367':[5,6,2,1,4,7,3,0],
          '0347v':[2,3,0,1,6,7,4,5],
          '0347':[6,7,3,2,5,4,0,1],
          '1256v':[0,1,2,3,4,5,6,7],
          '1256':[4,5,1,0,7,6,2,3]}

ca = vtk.vtkCellArray()
INEL = []
nPt = points.GetNumberOfPoints()
elist = []
parents = []
fsets = set()
for fg in femesh.face_groups:
    if fg.name in ['master','refine1']:
        if fg.name=='master':
            vert = True
        else:
            vert = False
        for kke in range(len(fg.element_faces)):
            face = fg.element_faces[kke]
            ke = face[0]-1
            elist.append(face[0])
            inel0 = femesh.vol.inel[ke]
            faceids = str()
            for kk,kn in enumerate(inel0):
                if kn in face[2]:
                    faceids += str(kk)
##            print(faceids)
            if len(faceids)==4:
        ##        print(faceids)
                if vert:
                    faceids += 'v'
                fsets.add(faceids)
                ie = iedict[faceids]
                inel = [inel0[kk] for kk in ie]
                # refined face: perpendicular to edge 1-2 on 2
                # refined direction: 2-3
                pts = [np.array([femesh.coords[kk][inel[kn]-1] for kk in range(3)]) for kn in range(8)]
                newpts = []
                # face0
                newpts.append(2/3*pts[1]+1/3*pts[5])
                newpts.append(1/3*pts[1]+2/3*pts[5])
                newpts.append(2/3*(1/3*pts[1]+2/3*pts[5])+1/3*(1/3*pts[0]+2/3*pts[4]))
                newpts.append(2/3*(2/3*pts[1]+1/3*pts[5])+1/3*(2/3*pts[0]+1/3*pts[4]))
                # face1
                newpts.append(2/3*pts[2]+1/3*pts[6])
                newpts.append(1/3*pts[2]+2/3*pts[6])
                newpts.append(2/3*(1/3*pts[2]+2/3*pts[6])+1/3*(1/3*pts[3]+2/3*pts[7]))
                newpts.append(2/3*(2/3*pts[2]+1/3*pts[6])+1/3*(2/3*pts[3]+1/3*pts[7]))
                newptids = []
                for kk in range(8):
                    id0 = pointLocator.IsInsertedPoint(newpts[kk])
                    if id0>-1:
                        newptids.append(id0)
                    else:
                        pointLocator.InsertPoint(pointLocator.GetPoints().GetNumberOfPoints(),newpts[kk])
                        newptids.append(pointLocator.GetPoints().GetNumberOfPoints()-1)
            ##            points.InsertNextPoint(newpts[kk])
                # back face
                hexa = vtk.vtkHexahedron()
                hexa.GetPointIds().SetId(0,inel[3]-1)
                hexa.GetPointIds().SetId(1,newptids[7])
                hexa.GetPointIds().SetId(2,newptids[6])
                hexa.GetPointIds().SetId(3,inel[7]-1)
                hexa.GetPointIds().SetId(4,inel[0]-1)
                hexa.GetPointIds().SetId(5,newptids[3])
                hexa.GetPointIds().SetId(6,newptids[2])
                hexa.GetPointIds().SetId(7,inel[4]-1)
                ca.InsertNextCell(hexa)
                parents.append(face[0])
                # lower face
                hexa = vtk.vtkHexahedron()
                hexa.GetPointIds().SetId(0,inel[3]-1)
                hexa.GetPointIds().SetId(1,inel[2]-1)
                hexa.GetPointIds().SetId(2,newptids[4])
                hexa.GetPointIds().SetId(3,newptids[7])
                hexa.GetPointIds().SetId(4,inel[0]-1)
                hexa.GetPointIds().SetId(5,inel[1]-1)
                hexa.GetPointIds().SetId(6,newptids[0])
                hexa.GetPointIds().SetId(7,newptids[3])
                ca.InsertNextCell(hexa)
                parents.append(face[0])
                # top face
                hexa = vtk.vtkHexahedron()
                hexa.GetPointIds().SetId(0,newptids[1])
                hexa.GetPointIds().SetId(1,newptids[2])
                hexa.GetPointIds().SetId(2,inel[4]-1)
                hexa.GetPointIds().SetId(3,inel[5]-1)
                hexa.GetPointIds().SetId(4,newptids[5])
                hexa.GetPointIds().SetId(5,newptids[6])
                hexa.GetPointIds().SetId(6,inel[7]-1)
                hexa.GetPointIds().SetId(7,inel[6]-1)
                ca.InsertNextCell(hexa)
                parents.append(face[0])
                # center piece
                hexa = vtk.vtkHexahedron()
                hexa.GetPointIds().SetId(0,newptids[0])
                hexa.GetPointIds().SetId(1,newptids[3])
                hexa.GetPointIds().SetId(2,newptids[2])
                hexa.GetPointIds().SetId(3,newptids[1])
                hexa.GetPointIds().SetId(4,newptids[4])
                hexa.GetPointIds().SetId(5,newptids[7])
                hexa.GetPointIds().SetId(6,newptids[6])
                hexa.GetPointIds().SetId(7,newptids[5])
                ca.InsertNextCell(hexa)
                parents.append(face[0])
mesh = vtk.vtkUnstructuredGrid()
mesh.SetPoints(points)
mesh.SetCells(vtk.VTK_HEXAHEDRON,ca)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('mesh.vtu')
writer.SetInputData(mesh)
writer.Write()

f = open(pathname+'/'+prob+'.inp')
eleprops = {}
for line in f:
    if '.i0g' in line:
        for ke in range(femesh.nVolumics):
            line = next(f)
            v = line.split(maxsplit=11)
            n = int(float(v[0]))
            if n in elist:
                eleprops[n] = v[-1]
        break
f.close()

f = open(pathname+'/M1224_3Dcomplet_empty.inp')
of = open(pathname+'/'+prob+'_ref.inp','w')

for line in f:
    if '.ing' in line:
        of.write(line)
        for kn in range(mesh.GetNumberOfPoints()):
            pt = mesh.GetPoint(kn)
            of.write('%i %1.10e %1.10e %1.10e 0\n'%(kn+1,pt[0],pt[1],pt[2]))
##        while len(line)>2:
##            line = next(f)
    elif '.i0g' in line:
        of.write(line)
        for ke in range(mesh.GetNumberOfCells()):
            cell = mesh.GetCell(ke)
            strinel = ''
            for kkn in range(8):
                strinel += '%i '%(cell.GetPointId(kkn)+1)
            of.write('%i %i B8 %s %s'%(ke+1,ke+1,strinel,eleprops[parents[ke]]))
##        while len(line)>2:
##            line = next(f)
    else:
        of.write(line)

f.close()
of.close()
                


