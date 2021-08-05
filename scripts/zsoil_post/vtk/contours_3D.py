# -*- coding: cp1252 -*-
import vtk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from zsoil_tools import vtktools

from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw



vrange = (0,0.04)
lut = vtktools.get_lut(lut_type='maps',vrange=vrange)
mlut = vtktools.get_lut(lut_type='mat')

etapes = {#1:u'Etat initial',
##          2:u'Réalisation pm, pieux sécants + berlinoise des carneaux, préfondés + barrettes butonnantes',
##          2.5:u'Réalisation 1er rang butons carneaux + partie Est dalle de couverture, ensuite 1ère étape d’excavation carneaux (niv. ~33.3 NGF)',
##          3:u'Réalisation 2ème rang butons carneaux, ensuite 2ème étape d’excavation carneaux (niv. ~29.1 NGF)',
##          4:u'Réalisation reste dalle de couverture, exc. S1 1ère étape',
##          4.5:u'Exc. S1 2ème étape',
##          5:u'Remblayage FUP',
##          6:r'Creusement du tunnel L14 en une étape, taux de déconfinement l=7%',
##          7:u'Réalisation dalle S1, exc. S2 1ère étape',
##          8:u'Exc. S2 2ème étape',
##          9:r'Creusement du tunnel L16/17 en une étape, taux de déconfinement l=7%',
##          9.5:u'Réalisation dalle S2, exc. S3 1ère étape',
##          10:u'Butonnage niv. 17.9 NGF',
##          10.5:u'Exc. S3 2ème étape',
##          11:u'Butonnage niv. 16.2 NGF',
##          11.5:u'Réalisation dalle S3 en laissant des brèches, exc. ff 1ère étape',
##          12:u'Butonnage niv 10.0 NGF',
          13:u'Exc. ff 2ème étape'}
##          14:u'Réalisation radier',
##          15:r'Dépose butons niv. 10.0 NGF centraux, creusement du tunnel L15 en une étape, taux de déconfinement l=7%',
##          16:u'Remplissage des brèches S3',
##          17:u'Dépose butons niv. 16.2 et 17.9 NGF',
##          18:u'Démolition barrettes butonnantes, application charges élé. non porteurs + charges mobiles'}

pathname = '../../pv'
prob = 'M1147_3D_SDP_coarse_v3_8(phasage)'

tx = []
f = open(pathname+'/'+prob+'.pvd')
for line in f:
    if 'shell.vtu' in line:
        t = float(line[line.find('timestep=')+9:].split('"')[1])
        if t in etapes.keys():
            tx.append(t)
f.close()
print(pathname+'/'+prob+'.pvd loaded')

for kt,t in enumerate(tx):
    tstr = str(int(t)).rjust(3,'0')+'_'+str(int((t-int(t))*100)).rjust(2,'0')

    meshes = []

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(pathname+'/'+prob+'_'+tstr+'_vol.vtu')
    reader.Update()
    meshes.append(reader.GetOutput())
    
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(pathname+'/'+prob+'_'+tstr+'_shell.vtu')
    reader.Update()
    shells = reader.GetOutput()
    sh2vol = vtk.vtkUnstructuredGrid()
    ploc = vtk.vtkPointLocator()
    bounds = meshes[0].GetBounds()
    ploc.InitPointInsertion(vtk.vtkPoints(),bounds)
    thick = shells.GetCellData().GetArray('THICK')
    mat = shells.GetCellData().GetArray('mat')
    EF = shells.GetCellData().GetArray('EF')
    for kc in range(shells.GetNumberOfCells()):
        pts0 = shells.GetCell(kc).GetPoints()
        pts01 = []
        for kk in range(4):
            pts01.append(np.array(pts0.GetPoint(kk)))
        normal = np.cross(pts01[1]-pts01[0],pts01[1]-pts01[3])
        normal /= np.linalg.norm(normal)
        th = thick.GetTuple1(kc)
        hexa = vtk.vtkHexahedron()
        for kv in range(4):
            apt = pts01[kv]-normal*th*0.5
            ptId = vtk.mutable(0)
            ploc.InsertUniquePoint(apt,ptId)
            hexa.GetPointIds().SetId(kv,ptId)
        for kv in range(4):
            apt = pts01[kv]+normal*th*0.5
            ptId = vtk.mutable(0)
            ploc.InsertUniquePoint(apt,ptId)
            hexa.GetPointIds().SetId(kv+4,ptId)
        sh2vol.InsertNextCell(hexa.GetCellType(),hexa.GetPointIds())
    sh2vol.SetPoints(ploc.GetPoints())
    sh2vol.GetCellData().AddArray(mat)
    sh2vol.GetCellData().AddArray(EF)
    sh2vol.GetCellData().AddArray(thick)
    
##    gf = vtk.vtkGeometryFilter()
##    gf.SetInputData(sh2vol)
####    gf.Update()
##    cpd = vtk.vtkCleanPolyData()
##    cpd.SetInputConnection(gf.GetOutputPort())
##    cpd.SetTolerance(0.01)
##    cpd.Update()
##    af = vtk.vtkAppendFilter()
##    af.AddInputConnection(cpd.GetOutputPort())
##    sh2vol2 = cpd.GetOutput()
##    print(sh2vol2.GetNumberOfCells())

##if True:
##    writer = vtk.vtkXMLUnstructuredGridWriter()
##    writer.SetFileName('vol.vtu')
##    writer.SetInputData(sh2vol2)
##    writer.Write()
        
    meshes.append(sh2vol)
    
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(pathname+'/'+prob+'_'+tstr+'_beam.vtu')
    reader.Update()
    meshes.append(reader.GetOutput())
    
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(pathname+'/'+prob+'_'+tstr+'_truss.vtu')
    reader.Update()
    meshes.append(reader.GetOutput())

    renderer = vtk.vtkRenderer()
    
    box = vtk.vtkBox()
    box.SetXMin(50.,-100.,-37.8)
    box.SetXMax(250.,100.,150.)
    plane = vtk.vtkPlane()
    plane.SetNormal(-1,0,0)
    plane.SetOrigin(162,0,0)
    
    for km,mesh in enumerate(meshes):
        if km==0:
            mesh.GetPointData().SetActiveScalars('DISP_TRA')
        else:
            mesh.GetCellData().SetActiveScalars('mat')

        cut = vtk.vtkExtractGeometry()
        cut.SetImplicitFunction(plane)
        cut.SetInputData(mesh)
        cut.SetExtractInside(False)
        cut.Update()
        extract = vtk.vtkExtractGeometry()
        extract.SetImplicitFunction(box)
        extract.SetExtractInside(False)
        extract.SetInputConnection(cut.GetOutputPort())
        extract.Update()
        if km==0:
            tf = vtk.vtkTransform()
            tf.Translate(-0.5,-0.2,-0.5)
            tff = vtk.vtkTransformFilter()
            tff.SetInputConnection(extract.GetOutputPort())
            tff.SetTransform(tf)
            tff.Update()
            extract = tff

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(extract.GetOutput())
        if km>0:
            mapper.SelectColorArray('mat')
            mapper.SetScalarRange((1,20))
            mapper.SetLookupTable(mlut)
        else:
            mapper.SelectColorArray('DISP_TRA')
            lut.SetRange(vrange)
            lut.SetVectorModeToMagnitude()
            mapper.SetScalarRange(vrange)
            mapper.SetLookupTable(lut)
            mapper.InterpolateScalarsBeforeMappingOn()
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        renderer.AddActor(actor)
        if km>1:
            actor.GetProperty().SetLineWidth(5)

##    ta = vtk.vtkTextActor()
##    ta.SetTextScaleModeToProp()
##    ta.SetDisplayPosition(10,10)
##    try:
##        ta.SetInput('T=%1.2f: %s'%(t,etapes[t]))
##    except:
##        ta.SetInput('T=%1.2f'%(t))
##    ta.SetMinimumSize(100,100)
##    tprop = ta.GetTextProperty()
##    tprop.SetFontSize(16)
##    tprop.SetFontFamilyToArial()
##    tprop.SetColor(0,0,0)
##    renderer.AddActor2D(ta)
    
    renderWin = vtk.vtkRenderWindow()
    renderWin.SetOffScreenRendering(1)
    renderWin.AddRenderer(renderer)
    renderWin.SetSize(1000,600)
    ##renderWin.SetAAFrames(5)
    renderWin.SetAlphaBitPlanes(1)

    center = meshes[0].GetCenter()
    camera=renderer.GetActiveCamera()
    camera.SetFocalPoint(center[0]+20,center[1],center[2])
    ####camera.SetClippingRange(100,10000)
    camera.SetPosition(center[0]+170,100,center[2]+200)
##    camera.SetViewAngle(10)
##    camera.SetDistance(10000)
    ##camera.SetViewUp(0,1,0)
##    camera.Zoom(0.1)
##    renderer.ResetCamera()
    renderer.SetBackground(1,1,1)
    renderWin.Render()

    screenshot=vtk.vtkWindowToImageFilter()
    screenshot.SetInput(renderWin)
    screenshot.SetMagnification(4)
    screenshot.SetInputBufferTypeToRGBA()
    screenshot.ReadFrontBufferOff()
    screenshot.Update()
    writer=vtk.vtkPNGWriter()
    writer.SetFileName(prob+'_maps_'+tstr+'.png')
    writer.SetInputData(screenshot.GetOutput())
    writer.Write()


    dpi = 96
    fig = plt.figure(figsize=(4,2100./dpi),dpi=dpi)
    bbox_props = dict(boxstyle='round',fc="w", ec="k", lw=2)

##    from matplotlib.patches import Rectangle
##    ax = fig.add_subplot(111)
##    nc = lut.GetNumberOfValues()
##    nval = nc
##    sp = 0.5
##    dh = min(0.8/(20*(1-sp)+(20+1)*sp)*(1-sp),1./(nval*(1-sp)+(nval+1)*sp)*(1-sp))
##    ds = dh/(1-sp)*sp
##    for kc in range(nval):
##        p = Rectangle((2*ds,1-(kc+1)*(dh+ds)),2-4*ds,dh,fc=lut.GetTableValue(kc),ec='k')
##        ax.add_patch(p)
##        ax.text(0.1,1-(kc+1)*(dh+ds),'%1.3g'%(float(kc)/nval),
##                va='center',size=20)
##    ax.annotate('Model created\nand computed\nwith ZSoil',size=24,
##                xy=(0.5,0.05),xycoords='axes fraction',ha='center',va='bottom')
##    ax.set_xlim(0,2)
##    ax.set_ylim(0,1)
##    ax.axis('off')
##    fig.tight_layout()
##    fig.savefig('legend.png')
##    leg = Image.open('legend.png')

    img = Image.open(prob+'_maps_'+tstr+'.png')
    img = img.crop((0,300,img.size[0]-200,img.size[1]))
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype(r'C:\Users\mprei\AppData\Local\Google\Chrome\User Data\Default\Extensions\eofcbnmajmjmplflapaojjnihcjkigck\12.0.911_0\common\ui\fonts\asp-open-sans\OpenSans-Regular.ttf',60)
    try:
        draw.text((100,img.size[1]-100),'T=%1.2f: %s'%(t,etapes[t]),font=font,fill=(0,0,0))
    except:
        draw.text((100,img.size[1]-100),'T=%1.2f'%(t),font=font,fill=(0,0,0))


    from zsoil_tools import postpro_lib as pl
    leg = pl.get_legend(lut,hfrac=0.8,vpad=0.1,hwratio=5,
                        label='Disp [m]')
    
    leg = leg.resize((int(float(leg.size[0])/leg.size[1]*img.size[1]),img.size[1]))
    img1 = Image.new('RGB',(img.size[0]+leg.size[0],img.size[1]))
    img1.paste(img,(0,0))
    img1.paste(leg,(img.size[0],0))
    
    img.save('disp_'+tstr+'.png')

    print('Phasage at T=%1.2f printed'%(t))

