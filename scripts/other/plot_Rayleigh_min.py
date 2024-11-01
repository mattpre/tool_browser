# -*- coding: cp1252 -*-
import numpy as np
import os,math
from numpy import linalg as la
from numpy import fft
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle



fqlogmin = -1.
fqlogmax = 2.
N = 100
T = [1./(10**(fqlogmin+(fqlogmax-fqlogmin)/(N-1)*kfq)) for kfq in range(N)]

# compute coefficients for x % damping at y Hz:
x = 1/1.2
y1 = 0.03
y2 = 0.005

ALPHA = [2*np.pi*y1/x,2*np.pi*y2/x]
BETA = [x*y1/2/np.pi,x*y2/2/np.pi]
labstr = [r'$\xi_{min}$=3% '+'@ %1.1f Hz:'%(1/x),
          r'$\xi_{min}$=0.5% '+'@ %1.1f Hz:'%(1/x)]

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

for kp in range(len(ALPHA)):
    xi = [100*(0.25/math.pi*ALPHA[kp]*t + math.pi*BETA[kp]/t) for t in T]
    ax.semilogx(T,xi,
                label=labstr[kp]+r' $\alpha=%1.4f$, $\beta=%1.2e$'%(ALPHA[kp],BETA[kp]))

##ax.add_patch(Rectangle((0.15,0),0.5-0.15,100,fc='grey',ec='none',alpha=0.5))
##ax.text(0.14,17,'Plateau BGK E',size=10,rotation=0)
##ax.add_patch(Rectangle((0.2,0),0.8-0.2,100,fc='grey',ec='none',alpha=0.5))
##ax.text(0.2,17,'Plateau BGK D',size=10,rotation=0)
##ax.add_patch(Rectangle((0.2,0),0.6-0.2,100,fc='grey',ec='none',alpha=0.5))
##ax.text(0.2,17,'Plateau classe de sol C',size=10,rotation=0)

ax.set_xlabel(u'Période [s]')
ax.set_title(u'Coefficient d\'amortissement hystérétique '+r'$\xi$ [%]')
ax.set_ylim([0,30])
ax.legend(prop={'size':'small'},loc=2)
ax.grid(which='both',axis='x')
ax.grid(which='major',axis='y')

fig.tight_layout()
fig.savefig('Rayleigh_min')

