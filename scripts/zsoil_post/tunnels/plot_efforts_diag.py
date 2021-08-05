import vtk
import numpy as np
import math
import matplotlib.pyplot as plt

from zsoil_tools import zsoil_results as zr

pathname = '..'
prob = 'M843_sortie3_mars18_2D_v1_2(sref,EI0)_sansCoin_surexcav'
##prob = 'M843_sortie3_mars18_2D_v1_1(sref)_l=70_E100'
##prob = 'M843_sortie3_mars18_2D_v1_1(sref)_l=70_E200'
##prob = 'M843_sortie3_mars18_2D_v1_1(1step)_25cm1step'

res = zr(pathname,prob)
res.read_rcf()
res.read_his()
tvect = [1.6,1.9,2]
tsteps = []
for kt,step in enumerate(res.steps):
    if step.conv_status==-1:
        for tx in tvect:
            if abs(step.time-tx)<1e-6:
                tsteps.append(kt)
                break
res.out_steps = tsteps
res.read_dat()
res.read_LTF()
res.read_s04()

cpos0 = []
for ke in range(res.nBeams):
    inel = res.beam.inel[ke]
    x = [res.coords[0][kn-1] for kn in inel]
    y = [res.coords[1][kn-1] for kn in inel]
    xm = 0.5*sum(x)
    ym = 0.5*sum(y)
    if np.angle(x[0]+y[0]*1j)>np.angle(x[1]+y[1]*1j):
        cpos0.append((xm+1000*ym,np.angle(xm+ym*1j),
                      (x,y)))
    else:
        cpos0.append((xm+1000*ym,np.angle(xm+ym*1j),
                      ([c for c in [x[1],x[0]]],
                       [c for c in [y[1],y[0]]])))
cpos0.sort(key=lambda v:v[1])
cpos = [v[0] for v in cpos0]
crds = [v[2] for v in cpos0]

M = [[[0,0] for kk in cpos] for k in range(3)]
N = [[0 for kk in cpos] for k in range(3)]
T = [[0 for kk in cpos] for k in range(3)]
for kkt,kt in enumerate(res.out_steps):
    step = res.steps[kt]
    for ke in range(res.nBeams):
        inel = res.beam.inel[ke]
        x = [res.coords[0][kn-1] for kn in inel]
        y = [res.coords[1][kn-1] for kn in inel]
        xm = 0.5*sum(x)
        ym = 0.5*sum(y)
        dl = ((x[1]-x[0])**2+(y[1]-y[0])**2)**0.5
        posind = cpos.index(xm+1000*ym)
        M[kkt][posind][0] = step.beam.moment[0][ke]+0.5*dl*step.beam.force[1][ke]
        M[kkt][posind][1] = step.beam.moment[0][ke]-0.5*dl*step.beam.force[1][ke]
        N[kkt][posind] = step.beam.force[0][ke]
        T[kkt][posind] = step.beam.force[1][ke]

sc = [5e-3,5e-4,2e-3]
AX = []
for kkt in range(len(res.out_steps)):
    step = res.steps[res.out_steps[kkt]]
    fig = plt.figure(figsize=(18,10))
    ax = plt.subplot(131)
    ax.set_title('Moment dans le cintre [kNm]')
    AX.append(ax)

    for ke in range(len(cpos)):
        p0 = np.array([crds[ke][0][0],crds[ke][1][0]])
        p1 = np.array([crds[ke][0][1],crds[ke][1][1]])
        v0 = (p1-p0)
        n = np.array([v0[1],-v0[0]])/np.linalg.norm(v0)
        ax.plot(crds[ke][0],crds[ke][1],'k')
        v11 = p0+n*sc[0]*M[kkt][ke][0]
        v12 = p1+n*sc[0]*M[kkt][ke][1]
        ax.plot([p0[0],v11[0],v12[0],p1[0]],
                [p0[1],v11[1],v12[1],p1[1]],'k')
        if 0.5*sum(M[kkt][ke])>0:
            val = max(M[kkt][ke])
        else:
            val = min(M[kkt][ke])
        ax.annotate('%1.0f'%(val),xy=0.5*(p0+p1)-0.4*n,
                    rotation=180+np.angle(n[0]+n[1]*1j,deg=1),
                    ha='center',va='center',size=10)

    ax = plt.subplot(132)
    ax.set_title('Effort normal dans le cintre [kN]')
    fig.text(0.01,0.01,prob,size=8)
    AX.append(ax)

    for ke in range(len(cpos)):
        p0 = np.array([crds[ke][0][0],crds[ke][1][0]])
        p1 = np.array([crds[ke][0][1],crds[ke][1][1]])
        v0 = (p1-p0)
        n = np.array([v0[1],-v0[0]])/np.linalg.norm(v0)
        ax.plot(crds[ke][0],crds[ke][1],'k')
        v11 = p0+n*sc[1]*N[kkt][ke]
        v12 = p1+n*sc[1]*N[kkt][ke]
        ax.plot([p0[0],v11[0],v12[0],p1[0]],
                [p0[1],v11[1],v12[1],p1[1]],'k')
        ax.annotate('%1.0f'%(N[kkt][ke]),xy=0.5*(p0+p1)-0.4*n,
                    rotation=180+np.angle(n[0]+n[1]*1j,deg=1),
                    ha='center',va='center',size=10)

    ax = plt.subplot(133)
    ax.set_title('Effort tranchant dans le cintre [kN]')
    AX.append(ax)

    for ke in range(len(cpos)):
        p0 = np.array([crds[ke][0][0],crds[ke][1][0]])
        p1 = np.array([crds[ke][0][1],crds[ke][1][1]])
        v0 = (p1-p0)
        n = np.array([v0[1],-v0[0]])/np.linalg.norm(v0)
        ax.plot(crds[ke][0],crds[ke][1],'k')
        v11 = p0+n*sc[2]*T[kkt][ke]
        v12 = p1+n*sc[2]*T[kkt][ke]
        ax.plot([p0[0],v11[0],v12[0],p1[0]],
                [p0[1],v11[1],v12[1],p1[1]],'k')
        ax.annotate('%1.0f'%(T[kkt][ke]),xy=0.5*(p0+p1)-0.4*n,
                    rotation=180+np.angle(n[0]+n[1]*1j,deg=1),
                    ha='center',va='center',size=10)

    for ax in AX:
        ax.grid('on')
        ax.set_aspect('equal')
        ax.set_xlim([0,4])
        ax.set_ylim([-1.5,5.5])
    
    fig.tight_layout()
    fig.savefig(prob+'_t=%i_%i'%(step.time,(step.time*10)%10))
    plt.close(fig)
