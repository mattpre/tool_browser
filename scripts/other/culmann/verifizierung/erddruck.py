import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from compute_Ea import *

H = 6

phi = np.radians(25)
delta = phi*2/3*0
gamma = 20


A = np.array([0,H])
C = [np.array([7,H])]
q = [np.array([0,0])]

figname = ''

tau = np.pi/2-delta
vect = np.array([-1,np.tan(tau-phi)])
vect /= np.linalg.norm(vect)
Ka = np.tan(np.pi/4-0.5*phi)**2
K0 = 1-np.sin(phi)
##K0 = (1-np.sin(phi))*(1+np.sin(beta))/np.cos(beta)


sc = 0.02
fig = plt.figure(figsize=(10,8))

ax = fig.add_subplot(111)

ax.plot([0,0],[0,A[1]],'k--')
pt0 = A
for kpt in range(len(C)):
    ax.plot([pt0[0],C[kpt][0]],[pt0[1],C[kpt][1]],'k')
    if np.linalg.norm(q[kpt])>0.01:
        ax.add_patch(Polygon([pt0,pt0+sc*q[kpt],C[kpt]+sc*q[kpt],C[kpt]]))
        ax.annotate('$q_k$=%1.1f kPa'%(np.linalg.norm(q[kpt])),
                    xy=(np.mean([pt0,pt0+sc*q[kpt],C[kpt]+sc*q[kpt],C[kpt]],0)),
                    ha='center')
    pt0 = C[kpt]

##ax.annotate('$\phi$ = %1.0f$^\circ$'%(np.degrees(phi)),xy=(-11.8,6.6))
##ax.annotate('$\delta$ = %1.2f$\phi$'%(delta/phia),xy=(-11.8,6.3))

EA = [0]
eA = [0]
Z = [A[1]]
Ind = [0]
PP = []
BB = []

nz = 11
for kz in range(1,nz):
    z = A[1]-A[1]/(nz-1)*kz
    Ea,ind,P,B = compute_Ea(A,C,q,phi,delta,gamma,tau,z,10)
    EA.append(Ea)
    Ind.append(ind)
    PP.append(P)
    BB.append(B)
##    E0.append(Ea/Kaa*K0a)
##    EAt.append(Ea-2*ca*(z1-z)*Kaa**0.5)
    Z.append(z)
##    print(EA[-2:])
##    print(Z[-2:])
    eA.append((EA[-1]-EA[-2])/(Z[-2]-Z[-1]))
##    eAt.append(max(0,(EAt[-1]-EAt[-2])/(Z[-2]-Z[-1])))
##    e0.append((E0[-1]-E0[-2])/(Z[-2]-Z[-1]))
        
##    print(x)

x = 0
ax.plot(np.array(eA)*sc-x,Z,'b',alpha=1,lw=1)
sc2 = 0.05
for kz in range(1,len(Z)):
    ax.plot([Ind[kz-1]-x,Ind[kz]-x],
            [Z[kz-1],Z[kz]],'k')
##    if abs(Z[kz]-8.55)<0.002:
    if kz in [1]:
        ax.plot([p[0]*sc2-x for p in PP[kz]],
                [p[1]*sc2+Z[kz] for p in PP[kz]],'r.-')
        dx = x+C[-1][0]
        ax.plot([-x,-x+dx],[Z[kz],Z[kz]+dx*np.tan(phi)],'k',lw=1,alpha=0.5)
        eam = 0
        ind = 0
        for kk in range(len(PP[kz])):
            ea = np.dot(PP[kz][kk]-BB[kz][kk],vect)
            if ea>eam:
                ind = kk
                eam = ea
        ax.plot([-x+BB[kz][ind][0]*sc2,-x+PP[kz][ind][0]*sc2],
                [Z[kz]+BB[kz][ind][1]*sc2,Z[kz]+PP[kz][ind][1]*sc2],'r')
                
        for kpt in range(len(C)):
            ax.plot([-x,C[kpt][0]],[Z[kz],C[kpt][1]],'k--',lw=1,alpha=0.3)
            
##    ax.annotate('%1.1f'%(eAe[ind]),xy=(eAe[ind]*sc-x,Z[ind]),va='top',
##                ha='right',xytext=(-2,-2),textcoords='offset points')
##    ind = 2
##    ax.annotate('%1.1f'%(eAe[ind]),xy=(eAe[ind]*sc-x,Z[ind]),
##                ha='left',xytext=(2,2),textcoords='offset points')
##    for kk in range(10,len(eAe)-10):
##        if (eAe[kk-1]-eAe[kk])/(eAe[kk]-eAe[kk+1])<1/20:
##            ax.annotate('%1.1f'%(eAe[kk]),xy=(eAe[kk]*sc-x,Z[kk]),va='top',
##                        ha='right',xytext=(-2,-2),textcoords='offset points')
##        elif (eAe[kk-1]-eAe[kk])/(eAe[kk]-eAe[kk+1])>12:
##            ax.annotate('%1.1f'%(eAe[kk]),xy=(eAe[kk]*sc-x,Z[kk]),va='bottom',
##                        ha='left',xytext=(2,2),textcoords='offset points')


ax.set_aspect('equal')
##ax.set_yticks([ax.get_ylim()[0]+ax.get_ylim()[1],
ax.grid(which='both',axis='both')
##ax.set_xlim([-1,5])
##ax.set_ylim([5,7])
##ax.legend([h0,h1,h2,h3],[u'Aktiver Erddruck ohne Kohäsion',
##                         u'Aktiver Erddruck mit Kohäsion',
##                         u'Erdruhedruck-Äquivalent nach SIA261',
##                         u'Erhöhter aktiver Erddruck'])

fig.tight_layout()
fig.savefig('Ea'+figname)
##fig.savefig('Ea_x=%1.0f'%(x))
plt.close(fig)
