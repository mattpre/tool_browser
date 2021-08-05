# -*- coding: cp1252 -*-
import numpy,os,math
from numpy import linalg as la
from numpy import fft
import matplotlib.pyplot as plt

from zsoil_tools import zsoil_results as zr
import pyrotd

def movingaverage(values,window):
    weigths = numpy.repeat(1.0, window)/window
    #including valid will REQUIRE there to be enough datapoints.
    #for example, if you take out valid, it will start @ point one,
    #not having any prior points, so itll be 1+0+0 = 1 /3 = .3333
    smas = numpy.convolve(values, weigths, 'valid')
    return smas # as a numpy array

def get_SIA_spect():
    # spektrum SIA 260:
    TT = [10**(kk*0.03)*0.01 for kk in range(101)]
    agd = 1.6
    q = 1.
    S = 1.4
    gf = 1.2
    TB = 0.15  #SIA261
    TC = 0.5  #SIA261
    TD = 2.0  #SIA261
##    TB = 0.15  #VD S8
##    TC = 0.5  #VD S8
##    TD = 2.0  #VD S8
    kB = [T>TB for T in TT].index(True)
    kC = [T>TC for T in TT].index(True)
    kD = [T>TD for T in TT].index(True)
    SIA = [gf*agd*S*(0.67+(2.5/q-0.67)*T/TB) for T in TT[:kB]]  #SIA261
    SIA.extend([2.5*gf*agd*S/q for T in TT[kB:kC]])  #SIA261
    SIA.extend([2.5*gf*agd*S/q*TC/T for T in TT[kC:kD]])  #SIA261
    SIA.extend([2.5*gf*agd*S/q*TC*TD/T**2 for T in TT[kD:]])  #SIA261

    return (TT,SIA)

pathname = '..'
pblist = ['1Dcol_dh2']

SA0 = [[] for k in pblist]

fx0 = [kf*0.05 for kf in range(1000)]

SA = []
for kf,prob in enumerate(pblist):
    res = zr(pathname,prob)
    res.read_rcf()
    res.read_his()
    tsteps = range(len(res.steps))
##        tsteps = range(500,1000)
    res.out_steps = tsteps
    res.read_dat()
    res.read_s00('/v+',res_type='accelerations')

    nd0 = [0,100]

    t0 = 0

    tx = []
    Ab = []
    At = []
    for kkt,kt in enumerate(res.out_steps):
        t = res.steps[kt].time
        tx.append(t-t0)
        Ab.append(res.steps[kt].nodal.a_disp[0][nd0[0]])
        At.append(res.steps[kt].nodal.a_disp[0][nd0[1]])

    fqlist = []
    fqlog = []
    fqlogmin = -1.
    fqlogmax = 2.
    N = 100
    for kfq in range(N):
        fqlist.append(10**(fqlogmin+(fqlogmax-fqlogmin)/(N-1)*kfq))
        fqlog.append((fqlogmin+(fqlogmax-fqlogmin)/(N-1)*kfq))
    T = [(1./f) for f in fqlist]
    #damping:
    xi = 0.05
    dt = tx[1]-tx[0]
    Sb = pyrotd.calc_spec_accels(dt,Ab,fqlist,xi)
    SA.append(Sb)
    St = pyrotd.calc_spec_accels(dt,At,fqlist,xi)
    SA.append(St)

if 1==1:
    fig = plt.figure(figsize=(6.5,4.5))
    AX = []

    TT,SIA = get_SIA_spect()

    ax = fig.add_subplot(111)
    for ka,a in enumerate(SA):
        ax.semilogx([1/v[0] for v in a],[v[1] for v in a],label=nd0[ka])
##    ax.semilogx(T,SA[0],'g-',label='Kontrollpunkt')
    ##ax.semilogx(TT,SIA,'k-.',label='SIA 261 Z3b, BGK E')
    ##ax.set_xlim(0,30)
    ##ax.set_ylim(0,2)
    ax.legend(prop={'size':'small'})

    ##ax.set_xticklabels(['%f'%(10**val for val in ax.get_xticks()])
    ##ax.set_xticklabels(['0.01','0.032','0.10','0.32','1.0','3.16','10.0'])
    ax.set_xlabel('Periode [s]')
    ax.set_ylabel('Spektralbeschleunigung [$m/s^2$]')
    AX.append(ax)


    for ax in AX:
        ax.grid(b=True,which='both',axis='x')
        ax.grid(b=True,which='major',axis='y')

    fig.tight_layout()
    fig.savefig('RS_'+prob)





