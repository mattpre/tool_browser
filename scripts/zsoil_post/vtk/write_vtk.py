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
prob = 'M1115_2020_3D_60cm_v2'

res = zr(pathname,prob)
res.read_rcf()
res.read_his()
tx = [0,2,3,4,5,6,7,8,9,10,11,12,13]
tsteps = []
for kt,step in enumerate(res.steps):
    if step.conv_status in [-1]:
        if abs(step.time%1)<1e-3:
            tsteps.append(kt)
res.out_steps = tsteps
res.read_dat()
res.read_s00()
for lab in res.ele_group_labels:
    if lab=='VOLUMICS':
        res.read_s01()  # volumics
    if lab=='SHELLS':
        res.read_s02()  # shells
##    elif lab=='TRUSSES':
##        res.read_s03()  # trusses
    elif lab=='BEAMS':
        res.read_s04()    # beams
##    elif lab=='CONTACT':
##        res.read_s07()

pvdFile = vtktools.create_pvd(res)
vtktools.write_vtu(res,beams=True,verbose=False,pvdFile=pvdFile)
##vtktools.write_vtu(res,trusses=True,verbose=False)
vtktools.write_vtu(res,vol=True,verbose=False,pvdFile=pvdFile)
vtktools.write_vtu(res,shells=True,verbose=False,pvdFile=pvdFile)
vtktools.save_pvd(pvdFile)


