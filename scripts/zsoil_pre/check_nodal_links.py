# -*- coding: cp1252 -*-
# @description Plotting evolution of nodal displacements.
# @output displacement plot
# @author Matthias Preisig
# @date 2017/08/18
import os,math
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle

from zsoil_tools import zsoil_results as zr

pathname = '.'
prob = '03-05 aprem_red_nlshell'


res = zr(pathname,prob)
res.read_dat()

nl1 = []
for kn in range(res.nNodes):
    if abs(res.coords[1][kn]-4.5)<2.5:
        if abs(res.coords[2][kn]+2)<0.02:
            nl1.append(kn+1)

for knl in range(res.nNodalLinks):
    klist = res.nodallinks[knl][1]-1
    kn = res.nodallinks[knl][0]
    dofs = res.nodallinks[knl][2]
    elist = res.lists[klist][0]
    ns = 0
    nv = 0
    nb = 0
    nt = 0
    for ke in elist:
        if ke in res.num_shells:
            ns += 1
        if ke in res.num_beams:
            nb += 1
        if ke in res.num_volumics:
            nv += 1
        if ke in res.num_trusses:
            nt += 1
##    if nv>0 and 'LList' in res.lists[klist][1]:
##        print('%s: sh=%i, vol=%i, be=%i, tr=%i'%(res.lists[klist][1],ns,nv,nb,nt))
    if kn in nl1 and 'LList' in res.lists[klist][1]:
        print('%s: kn=%i, dofs=%i, sh=%i, vol=%i, be=%i, tr=%i, (%1.3f,%1.3f,%1.3f)'%(res.lists[klist][1],kn,dofs,ns,nv,nb,nt,res.coords[0][kn],res.coords[1][kn],res.coords[2][kn]))
            
        

            




