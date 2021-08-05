# @description Exporting zsoil results to vtu
# @input zsoil results
# @output vtu unstructured grid
# @author Matthias Preisig
# @date 2017/10/10
import numpy,os
from numpy import linalg as la


from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools

pathname = '..'
pblist = ['M1292_3D_v1']

for prob in pblist:
    res = zr(pathname,prob)
    res.read_rcf()
    res.read_his()
    tx = [1.5,3,4]
    tsteps = []
    for kt,step in enumerate(res.steps):
        if step.conv_status in [-1]:
##            if step.time in tx:
            tsteps.append(kt)
    tsteps.append(len(res.steps)-1)
    res.out_steps = tsteps
    res.read_dat()
    res.read_s00()
    for lab in res.ele_group_labels:
        if lab=='VOLUMICS':
            res.read_s01()  # volumics
        elif lab=='SHELLS':
            res.read_s02()  # shells
        elif lab=='TRUSSES':
            res.read_s03()  # trusses
##        elif lab=='BEAMS':
##            res.read_s04()    # beams
        elif lab=='CONTACT':
            res.read_s07()

    step0 = res.steps[tsteps[0]]

    pvdFile = vtktools.create_pvd(res)
##    vtktools.write_vtu(res,beams=True,verbose=False,pvdFile=pvdFile)#,refstep=step0,tsteps=tsteps[1:])
    vtktools.write_vtu(res,trusses=True,verbose=False,pvdFile=pvdFile)#
    vtktools.write_vtu(res,vol=True,verbose=False,pvdFile=pvdFile)#,refstep=step0,tsteps=tsteps[1:])
    vtktools.write_vtu(res,shells=True,verbose=False,pvdFile=pvdFile)#,refstep=step0,tsteps=tsteps[1:])
##    vtktools.write_vtu(res,cnt=True,verbose=False,pvdFile=pvdFile)#,refstep=step0,tsteps=tsteps[1:])
    vtktools.save_pvd(pvdFile)


