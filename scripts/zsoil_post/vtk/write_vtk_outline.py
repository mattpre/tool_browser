# @description Exporting outline (boundary faces) of zsoil results to vtu
# @input zsoil results
# @output vtu unstructured grid
# @author Matthias Preisig
# @date 2017/10/10

import numpy as np

from zsoil_tools import zsoil_results as zr
from zsoil_tools import vtktools


pathname = r'\\192.168.1.51\Mandats sur H RAID0\M1010_Tourbillon\stab_panneau'
prob = 'M1010_stabPann_m2_renfLat'

res = zr(pathname,prob)
res.read_rcf()
res.read_his()
tx = [67]
tsteps = []
for kt,step in enumerate(res.steps):
    if step.conv_status in [-1]:
        if step.time in tx:
            tsteps.append(kt)
res.out_steps = tsteps
res.read_dat()
res.read_s00()
for lab in res.ele_group_labels:
    if lab=='VOLUMICS':
        res.read_s01()  # volumics
##    elif lab=='SHELLS':
##        res.read_s02()  # shells
##    elif lab=='TRUSSES':
##        res.read_s03()  # trusses
##    elif lab=='BEAMS':
##        res.read_s04()    # beams
##    elif lab=='CONTACT':
##        res.read_s07()


##vtktools.write_vtu(res,beams=True,verbose=False)
##vtktools.write_vtu(res,trusses=True,verbose=False)
vtktools.write_vtu(res,vol=True,verbose=False,outline=True)
##vtktools.write_vtu(res,shells=True,verbose=False)

