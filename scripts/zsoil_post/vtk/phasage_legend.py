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

matdict = {1:u'Remblais',
           3:u'MarnesàPholadomyes2',
           4:u'Sables verts 1',
           5:u'Calcaire SO2',
           6:u'Calcaire SO3',
           7:u'Sables de Beauchamp1',
           8:u'MarnesCaillasses3',
           9:u'MarnesCaillasses4',
           10:u'CalcaireGrossier',
           11:u'PM',
           12:u'voussoirs',
           13:u'pieux sécants',
           14:u'paroi berlinoise',
           15:u'dalles',
           16:u'colonnes',
           17:u'barrettes',
           18:u'radier',
           19:u'longrines',
           21:u'seepage',
           22:u'cntPM',
           25:u'qs',
           26:u'butons',
           27:u'semelles bâtiments',
           28:u'butons 1 carneaux',
           29:u'butons 2 carneaux'}

lut = vtktools.get_lut(lut_type='mat',maxind=29)

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
        mapper.SelectColorArray('mat')
        mapper.SetScalarRange((1,lut.GetNumberOfColors()))
        mapper.SetLookupTable(lut)
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
    writer.SetFileName(prob+'_phasage_'+tstr+'.png')
    writer.SetInputData(screenshot.GetOutput())
    writer.Write()

    ## dalles et etayage:
    fig = plt.figure(figsize=(11,8.5))
    bbox_props = dict(boxstyle='round',fc="w", ec="k", lw=2)

    ax = fig.add_subplot(111)
    
    nivs = [37.5,29.05,22.05,15.2,5.3]
    nstr = ['Dalle de couverture',
            'Dalle S1','Dalle S2','Dalle S3','Radier']
    nivt = [0,4,7,9.5,13]
    nind = max([n for n,i in enumerate(nivt) if i<t])

    trusses = meshes[3]
    trtype = [('b',1),('m',2),('r',3)]
    trnivlist = set()
    for kc in range(trusses.GetNumberOfCells()):
        cell = trusses.GetCell(kc)
        pts = [trusses.GetPoint(cell.GetPointId(kp)) for kp in range(2)]
        xm = 0.5*sum([pt[0] for pt in pts])
        ym = 0.5*sum([pt[1] for pt in pts])
        zm = 0.5*sum([pt[2] for pt in pts])
        if abs(ym-17.9)<1e-3:
            ax.plot([pt[0] for pt in pts],[-pt[2] for pt in pts],'b',linewidth=1)
            trnivlist.add(0)
        elif abs(ym-16.2)<1e-3:
            ax.plot([pt[0] for pt in pts],[-pt[2] for pt in pts],'r',linewidth=3)
            trnivlist.add(1)
        elif abs(ym-10.0)<1e-3:
            ax.plot([pt[0] for pt in pts],[-pt[2] for pt in pts],'g',linewidth=6)
            trnivlist.add(2)
    if 0 in list(trnivlist):
        ax.plot([40,45],[30,30],'b',linewidth=1)
        ax.annotate('Butons niv. 17.9 NGF',xy=(48,30),va='center',size=20)
    if 1 in list(trnivlist):
        ax.plot([40,45],[26,26],'r',linewidth=3)
        ax.annotate('Butons niv. 16.2 NGF',xy=(48,26),va='center',size=20)
    if 2 in list(trnivlist):
        ax.plot([40,45],[22,22],'g',linewidth=6)
        ax.annotate('Butons niv. 10.0 NGF',xy=(48,22),va='center',size=20)

    patches = []
    cvect = []
    for kc in range(shells.GetNumberOfCells()):
        cell = shells.GetCell(kc)
        pts = [shells.GetPoint(cell.GetPointId(kp)) for kp in range(4)]
        xm = 0.25*sum([pt[0] for pt in pts])
        ym = 0.25*sum([pt[1] for pt in pts])
        zm = 0.25*sum([pt[2] for pt in pts])
        if abs(ym-nivs[nind])<1e-3:
##            ax.plot([pt[0] for pt in pts],[-pt[2] for pt in pts],'k')
            xx = [pt[0] for pt in pts]
            zz = [-pt[2] for pt in pts]
            patches.append(Polygon(np.array([xx,zz]).T))
            cvect.append(thick.GetTuple1(kc))
    pc = PatchCollection(patches,edgecolors='none',facecolors='r',alpha=0.5)
##    pc.set_array(np.array(cvect))
    ax.add_collection(pc)

    ax.annotate(nstr[nind],xy=(0.5,0.5),xycoords='axes fraction',
                ha='center',va='center',size=28)

    ax.set_xlim([-2,125])
    ax.set_ylim([-2,82])
    ax.axis('off')
    ax.set_aspect('equal')

    fig.tight_layout()
    fig.savefig('dalle_tmp.png',transparent=True)
    plt.close(fig)

    img = Image.open(prob+'_phasage_'+tstr+'.png')
    img = img.crop((0,300,img.size[0]-400,img.size[1]))
    dalles = Image.open('dalle_tmp.png')
    img.paste(dalles,(int(0.29*img.size[0]),img.size[1]-dalles.size[1]),dalles)

    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype(r'C:\Users\mprei\AppData\Local\Google\Chrome\User Data\Default\Extensions\eofcbnmajmjmplflapaojjnihcjkigck\12.0.911_0\common\ui\fonts\asp-open-sans\OpenSans-Regular.ttf',60)
    try:
        draw.text((100,img.size[1]-100),'T=%1.2f: %s'%(t,etapes[t]),font=font,fill=(0,0,0))
    except:
        draw.text((100,img.size[1]-100),'T=%1.2f'%(t),font=font,fill=(0,0,0))


    from zsoil_tools import postpro_lib as pl
    nval = lut.GetNumberOfColors()
    leg = pl.get_legend(lut,categories=matdict,hfrac=0.8,vpad=0.1,hwratio=5,
                        label='Disp [m]')
    
    leg = leg.resize((int(float(leg.size[0])/leg.size[1]*img.size[1]),img.size[1]))
    img1 = Image.new('RGB',(img.size[0]+leg.size[0],img.size[1]))
    img1.paste(img,(0,0))
    img1.paste(leg,(img.size[0],0))
    
    img1.save('phasage_'+tstr+'.png')

    print('Phasage at T=%1.2f printed'%(t))

