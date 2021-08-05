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
          5:u'Tragwerkslasten'}

pathname = '../pv'
pblist = ['M1013_O2021_v1_pilesV1017_cent1_L1Ost',
          'M1013_O2021_v1_pilesV1017_cent1_L2Ost',
          'M1013_O2021_v1_piles_red1_L1Ost',
          'M1013_O2021_v1_noPiles_L1Ost',
          'M1013_O_v8_pilesV1017_cent2_l171110']

for prob in pblist:

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
        tstr = vtktools.get_tstr(t,t0=4)
    ##    tstr = str(int(t)).rjust(3,'0')+'_'+str(int((t-int(t))*100)).rjust(2,'0')

        meshes = []

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(pathname+'/'+prob+'_'+tstr+'_vol.vtu')
        reader.Update()
        meshes.append(reader.GetOutput())

        arr = meshes[0].GetPointData().GetArray('DISP_TRA')
##        vrange = (0,arr.GetMaxNorm())
        vrange = (arr.GetRange(1)[0],0)
        lut = vtktools.get_lut(lut_type='maps',vrange=vrange)
        mlut = vtktools.get_lut(lut_type='mat')

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
            
    ##    meshes.append(sh2vol)
    ##    
    ##    reader = vtk.vtkXMLUnstructuredGridReader()
    ##    reader.SetFileName(pathname+'/'+prob+'_'+tstr+'_beam.vtu')
    ##    reader.Update()
    ##    meshes.append(reader.GetOutput())
        
    ##    reader = vtk.vtkXMLUnstructuredGridReader()
    ##    reader.SetFileName(pathname+'/'+prob+'_'+tstr+'_truss.vtu')
    ##    reader.Update()
    ##    meshes.append(reader.GetOutput())

        renderer = vtk.vtkRenderer()
        
        box = vtk.vtkBox()
        box.SetXMin(200.,-100.,-330.0)
        box.SetXMax(2500.,400.,2500.)
        plane = vtk.vtkPlane()
        plane.SetNormal(-1,0,0)
        plane.SetOrigin(0,0,0)
        
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
    ##        extract = vtk.vtkExtractGeometry()
    ##        extract.SetImplicitFunction(box)
    ##        extract.SetExtractInside(False)
    ##        extract.SetInputConnection(cut.GetOutputPort())
    ##        extract.Update()
            extract = vtk.vtkClipDataSet()
            extract.SetClipFunction(box)
            extract.SetInsideOut(False)
##            extract.SetInputConnection(cut.GetOutputPort())
            extract.SetInputData(mesh)
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
##                lut.SetVectorModeToMagnitude()
                lut.SetVectorModeToComponent()
                lut.SetVectorComponent(1)
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
        renderer.RemoveAllLights()
        lightkit = vtk.vtkLightKit()
        lightkit.SetKeyToFillRatio(2)
        lightkit.AddLightsToRenderer(renderer)
        
        renderWin = vtk.vtkRenderWindow()
        renderWin.SetOffScreenRendering(1)
        renderWin.AddRenderer(renderer)
        renderWin.SetSize(800,500)
        ##renderWin.SetAAFrames(5)
        renderWin.SetAlphaBitPlanes(1)

        center = meshes[0].GetCenter()
        camera=renderer.GetActiveCamera()
        camera.SetFocalPoint(center[0],center[1]+80,center[2])
        ####camera.SetClippingRange(100,10000)
        camera.SetPosition(center[0]+120,520,center[2]+120)
    ##    camera.SetViewAngle(10)
    ##    camera.SetDistance(10000)
        ##camera.SetViewUp(0,1,0)
    ##    camera.Zoom(0.1)
    ##    renderer.ResetCamera()
        renderer.SetBackground(1,1,1)
        renderWin.Render()

        screenshot=vtk.vtkWindowToImageFilter()
        screenshot.SetInput(renderWin)
        screenshot.SetScale(4)
    ##    screenshot.SetMagnification(4)
        screenshot.SetInputBufferTypeToRGBA()
        screenshot.ReadFrontBufferOff()
        screenshot.Update()
        writer=vtk.vtkPNGWriter()
        writer.SetFileName(prob+'_maps_'+tstr+'.png')
        writer.SetInputData(screenshot.GetOutput())
        writer.Write()


##        dpi = 96
##        fig = plt.figure(figsize=(4,2100./dpi),dpi=dpi)
##        bbox_props = dict(boxstyle='round',fc="w", ec="k", lw=2)

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
    ##    img = img.crop((0,300,img.size[0]-200,img.size[1]))
        draw = ImageDraw.Draw(img)
        draw.rectangle((50,10,1000,200),fill='white')
##        draw.rectangle((10,img.size[1]-20,400,img.size[1]),fill='white')
        font = ImageFont.truetype("arial",60)
        fonts = ImageFont.truetype("arial",20)
        try:
            draw.text((100,50),'T=%1.2f: %s'%(t,etapes[t]),font=font,fill=(0,0,0))
        except:
            draw.text((100,50),'T=%1.2f'%(t),font=font,fill=(0,0,0))
        draw.text((10,img.size[1]-20),'%s'%(prob),font=fonts,fill='white')


        from zsoil_tools import postpro_lib as pl
        leg = pl.get_legend(lut,hfrac=1.0,vpad=0.1,hwratio=5,
                            label='Setzung [m]')
        
        leg = leg.resize((int(float(leg.size[0])/leg.size[1]*img.size[1]),img.size[1]))
        img1 = Image.new('RGB',(img.size[0]+leg.size[0],img.size[1]))
        img1.paste(img,(0,0))
        img1.paste(leg,(img.size[0],0))
        
        img1.save('disp_%s_%s.png'%(prob,tstr))

        print('Phasage at T=%1.2f printed'%(t))

