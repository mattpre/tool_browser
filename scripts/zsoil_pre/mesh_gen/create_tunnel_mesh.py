# @create layered mesh between two opposing faces. used for tunnel support with different layers (original concrete, demolition, shotcrete, sealing etc.)
# @date 2018/10/24
# @author Matthias Preisig

import numpy as np
import math

f = open('M1135_Var0.inp')

thick = 0.34

var = 'Var1'
rlist = [0.06,0.10,0.109,0.137,0.14,0.3,thick]
split = [8,6,2,6,1,8,2]
EFs = [2,3,6,6,6,10,7]
rm1 = [2,2,2,2,2,2,3]
rm2 = [0,4,6,4,5,3,0]
load_kc = [0,5]

##var = 'Var2'
##rlist = [0.09,0.13,0.139,0.142,0.3,thick]
##split = [6,4,2,1,8,2]
##EFs = [2,3,6,6,10,7]
##rm1 = [2,2,2,2,2,3]
##rm2 = [0,4,6,5,3,0]
##load_kc = [0,4]

of = open('M1135_%s_meshgen.inp'%(var),'w')
nilist = []
nelist = []

c0 = np.array([0,2.0703378,0])

for line in f:
    if '.ing' in line:
        of.write(line)
        line = f.next()
        while len(line)>1:
            of.write(line)
            v = line.split()
            pt = np.array([float(v[1]),float(v[2]),0])
            if abs(np.linalg.norm(pt-c0)-4.96)<1e-4 and pt[1]>-1e-6:
                nilist.append((int(v[0]),pt,np.angle((pt[0]-c0[0])*1j+(pt[1]-c0[1]))))
            elif abs(np.linalg.norm(pt-c0)-5.3)<1e-4 and pt[1]>-1e-6:
                nelist.append((int(v[0]),pt,np.angle((pt[0]-c0[0])*1j+(pt[1]-c0[1]))))
            line = f.next()
        last_nd = int(v[0])
        nilist.sort(key=lambda v:v[2])
        nelist.sort(key=lambda v:v[2])
        nelist = nelist[::2]
        knd = 0
        r0 = 0
        for kc in range(len(rlist)):
            for kk in range(split[kc]):
                if not (kc==len(rlist)-1 and kk==split[kc] and kc==0 and kk==0):
                    for kn in range(len(nilist)):
                        knd += 1
                        pte = nelist[kn][1]
                        pti = nilist[kn][1]
                        pt = pte+((r0+(rlist[kc]-r0)*float(kk+1)/split[kc])/thick)*(pti-pte)
                        of.write('%i %1.12e %1.12e %1.12e 0\n'%(last_nd+knd,pt[0],pt[1],0))
            r0 = rlist[kc]
        of.write(line)
    elif '.i0g' in line:
        of.write(line)
        line = f.next()
        while len(line)>1:
            of.write(line)
            last_ele = int(line.split()[0])
            last0 = last_ele
            line = f.next()
        for ke in range(len(nilist)-1):
            ne0 = nelist[ke][0]
            ne1 = nelist[ke+1][0]
            nc = 0
            for kc in range(len(rlist)):
                ef = EFs[kc]
                RM1 = rm1[kc]
                RM2 = rm2[kc]
                if kc==1:
                    if 0.5*(nelist[ke][2]+nelist[ke+1][2])<-1.1595 or 0.5*(nelist[ke][2]+nelist[ke+1][2])>1.1595:
                        ef = EFs[0]
                        RM1 = rm1[0]
                        RM2 = rm2[0]
                for kk in range(split[kc]):
                    if not (kc==len(rlist)-1 and kk==split[kc]-1):
                        ni0 = last_nd+1+nc*(len(nelist))+ke
                        ni1 = last_nd+2+nc*(len(nelist))+ke
                    else:
                        ni0 = nilist[ke][0]
                        ni1 = nilist[ke+1][0]
                    of.write('%i %i Q4 %i %i %i %i %i %i %i %i %i %i %i 0\n'%\
                             (last_ele+1,last_ele+1,ne1,ne0,ni0,ni1,1,1,1,RM1,RM2,ef,0))
                    last_ele += 1
                    ne0 = ni0
                    ne1 = ni1
                    nc += 1
    elif '.isl' in line:
        of.write(line)
        kl = 0
        last_ele = last0
        for ke in range(len(nilist)-1):
            for kc in range(len(rlist)):
                if kc==load_kc[0] and (nelist[ke][1][0]>-2.73 and nelist[ke][1][0]<4.79):
                    of.write('%i %i %i %i %i 0 0\n'%(kl+1,last_ele+1,1,2,2))
                    of.write('unif_1\n')
                    of.write(' %i %1.12e %1.12e %1.12e\n'%(0,0,0,0))
                    of.write(' %i %1.12e %1.12e %1.12e\n'%(1,-31,0,0))
                    kl += 1
                elif kc==load_kc[1] and (nelist[ke][1][0]>-2.82 and nelist[ke][1][0]<4.87):
                    of.write('%i %i %i %i %i 0 0\n'%(kl+1,last_ele+1,1,5,2))
                    of.write('unif_2\n')
                    of.write(' %i %1.12e %1.12e %1.12e\n'%(0,0,0,0))
                    of.write(' %i %1.12e %1.12e %1.12e\n'%(1,-31,0,0))
                    kl += 1
                for kk in range(split[kc]):
                    last_ele += 1
    elif '.gsl' in line:
        of.write(line)
        kl = 0
        last_ele = last0
        for ke in range(len(nilist)-1):
            for kc in range(len(rlist)):
                if kc==load_kc[0] and (nelist[ke][1][0]>-2.73 and nelist[ke][1][0]<4.79):
                    of.write('%i UNI_LOAD 0 %i 0 %i 0\n'%(kl+1,2,1))
                    of.write('unif_1\n')
                    of.write('%1.12e %1.12e %1.12e\n'%(0,-31,0))
                    of.write('%i 1\n'%(last_ele+1))
                    kl += 1
                elif kc==load_kc[1] and (nelist[ke][1][0]>-2.82 and nelist[ke][1][0]<4.87):
                    of.write('%i UNI_LOAD 0 %i 0 %i 0\n'%(kl+1,5,1))
                    of.write('unif_2\n')
                    of.write('%1.12e %1.12e %1.12e\n'%(0,-31,0))
                    of.write('%i 1\n'%(last_ele+1))
                    kl += 1
                for kk in range(split[kc]):
                    last_ele += 1
    else:
        of.write(line)

of.close()
f.close()
