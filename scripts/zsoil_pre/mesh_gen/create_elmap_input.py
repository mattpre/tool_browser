# @description deformes regular mesh to circular tunnel
# @date 2017/1/26
# @author Matthias Preisig

import cmath,math


def on_edge(c,edge):
    tol = 1.e-6
    a = (edge[0][0],edge[1][0])
    b = (edge[0][1],edge[1][1])
    
    crossproduct = (c[1]-a[1])*(b[0]-a[0])-(c[0]-a[0])*(b[1]-a[1])
    if abs(crossproduct)>tol:
        return False,0

    dotproduct = (c[0]-a[0])*(b[0]-a[0])+(c[1]-a[1])*(b[1]-a[1])
    if dotproduct<0:
        return False,0

    squaredlengthba = (b[0]-a[0])*(b[0]-a[0])+(b[1]-a[1])*(b[1]-a[1])
    if dotproduct>squaredlengthba:
        return False,0

    return True,dotproduct/squaredlengthba

crd = [[],[]]

f = open('extrusion_tunnel_fine_template.inp')
for line in f:
    if '.ing' in line:
        line = f.next()
        while not line[0]=='.' and len(line)>1:
            v = line.split()
            for kk in range(2):
                crd[kk].append(float(v[kk+1]))
            line = f.next()
f.close()

edgenodes = [[1.75,1.75,3.25,3.25],
             [3.25,1.75,1.75,3.25]]
edges = [[0,1],[1,2],[2,3],[3,0]]

# find nodes on edges:
en = []
for kn in range(len(crd[0])):
    for ke in range(len(edges)):
        e = [[edgenodes[0][edges[ke][0]],edgenodes[0][edges[ke][1]]],
             [edgenodes[1][edges[ke][0]],edgenodes[1][edges[ke][1]]]]
        res = on_edge((crd[0][kn],crd[1][kn]),e)
        if res[0]:
            en.append((kn,ke,res[1]))
            break

c = (2.5,2.5)
r = 0.9
dx = []
for ke in range(len(en)):
    crd0 = (crd[0][en[ke][0]],crd[1][en[ke][0]])
    # project on circle with center c and radius r
    ptc = (crd0[0]-c[0])+(crd0[1]-c[1])*1j
    ang = cmath.phase(ptc)
    crd1 = (c[0]+r*cmath.cos(ang).real,c[1]+r*cmath.sin(ang).real)
    dx.append((crd1[0]-crd0[0],crd1[1]-crd0[1]))


f = open('extrusion_tunnel_fine_template.inp')
of = open('extrusion_tunnel_def.inp','w')
of.write(f.next())
line = f.next()
v = line.split()
string = ''
for kk in range(len(v)):
    if kk==6:
        string += '%i '%(int(v[kk])+len(en))
    else:
        string += '%s '%(v[kk])
of.write(string+'\n')
for line in f:
    if '.inb' in line:
        of.write(line)
        count = 0
        while len(line)>1:
            line = f.next()
            if len(line)>1:
                count += 1
                of.write(line)
        for kn in range(len(en)):
            of.write('%i %i 1 1 %1.8e 0 1 0 1 %1.8e 0 1 0 0 %1.8e 0 0 0 0\n'%
                     (kn+1+count,en[kn][0]+1,dx[kn][0],dx[kn][1],0))
        of.write('\n')
    elif '.inc' in line:
        of.write(line)
        count = 0
        while len(line)>1:
            for kk in range(21):
                line = f.next()
                if len(line)<3:
                    break
                of.write(line)
            if len(line)>2:
                count += 1
        for kn in range(len(en)):
            of.write('%i %i 0 0 0\n'%(kn+1+count,en[kn][0]+1))
            of.write('1 1 0 0 0 0\n')
            of.write('1 %1.8e 0 1 0\n'%(dx[kn][0]))
            of.write('0 %1.8e 0 0\n'%(0))
            of.write('0 %1.8e 0 0\n'%(0))
            of.write('1 %1.8e 0 1 0\n'%(dx[kn][1]))
            of.write('0 %1.8e 0 0\n'%(0))
            of.write('0 %1.8e 0 0\n'%(0))
            for kk in range(4):
                of.write('0 %1.8e 0 0 0\n'%(0))
                of.write('0 %1.8e 0 0\n'%(0))
                of.write('0 %1.8e 0 0\n'%(0))
            of.write('No name\n')
        of.write('\n')
    else:
        of.write(line)
f.close()
of.close()

##import os
##cd = os.getcwd()
##zspath = '\"C:\\Program Files\\ZSoil\\ZSoil 2016 v16.03 x64'
##os.system(zspath+'/Z_Soil.exe\" '+cd+'\\extrusion_tunnel_def.inp /E')



import sys
sys.path.append('D:/Mandats/python/zsoil_tools')
import zsoil_results as zr

pathname = '.'
prob = 'extrusion_tunnel_def'

res = zr.zsoil_results(pathname,prob)
res.read_rcf()
res.read_his()

tsteps = []
for kt,step in enumerate(res.steps):
    if step.time==1 and step.sf==0 and step.conv_status==-1:
        tsteps.append(kt)

res.out_steps = tsteps
res.read_dat()
res.read_s00()




f = open('extrusion_tunnel_def.inp')
of = open('extrusion_tunnel_def_mod.inp','w')

for line in f:
    if '.ing' in line:
        of.write(line)
        disp = res.steps[1].nodal.disp
        for kn in range(res.nNodes):
            line = f.next()
            of.write('%i %1.8e %1.8e %1.8e 0\n'%
                     (kn+1,res.coords[0][kn]+disp[0][kn],
                      res.coords[1][kn]+disp[1][kn],0))
        of.write('\n')
    else:
        of.write(line)
f.close()
of.close()
