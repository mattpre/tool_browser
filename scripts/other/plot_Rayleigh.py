# -*- coding: cp1252 -*-
import numpy,os,math
from numpy import linalg as la
from numpy import fft
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle



fqlogmin = -1.
fqlogmax = 2.
N = 100
T = [1./(10**(fqlogmin+(fqlogmax-fqlogmin)/(N-1)*kfq)) for kfq in range(N)]

ALPHA = [0.05,1.176,0.05,0.232]
BETA = [0,0.001592,0.001592,0.008]
labstr = ['Boden Geomod',
          'Boden HSR',
          r'Boden HSR $\beta$',
          u'Stützwand HSR']
ALPHA = [0.966643893412244]
BETA = [0.00183640318952187]
labstr = ['Amortissement de Rayleigh']

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

for kp in range(len(ALPHA)):
    xi = [100*(0.25/math.pi*ALPHA[kp]*t + math.pi*BETA[kp]/t) for t in T]
    ax.semilogx(T,xi,
                label=labstr[kp]+r' $\alpha=%1.3f$, $\beta=%1.1e$'%(ALPHA[kp],BETA[kp]))

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
ax.grid(b=True,which='both',axis='x')
ax.grid(b=True,which='major',axis='y')

fig.tight_layout()
fig.savefig('Rayleigh')

