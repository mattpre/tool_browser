# -*- coding: cp1252 -*-
import numpy as np
import os,math
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle





fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

z = np.linspace(-10,0,100)
p = [v*10+80 for v in z]

for alpha in [2,5,10,45]:
    Sr = 0.03
    gf = 10.
    S = [1 if v<=0 else Sr+(1-Sr)/(1+(alpha*v/gf)**2)**0.5 for v in p]

    ax.plot(S,z,label='$\\alpha=$%1.1f $m^{-1}$'%(alpha))
ax.plot([0,1],[-8,-8],'b--')
ax.annotate('Water table',xy=(0.5,-8),ha='center',va='top')
ax.plot([Sr,Sr],[min(z),max(z)],'k--')
ax.annotate('Sr=%1.2f'%(Sr),xy=(Sr,-5),rotation=90,ha='right',va='center')

ax.set_xlabel(u'Saturation ratio [-]')
ax.set_ylabel(u'Depth [m]')
ax.legend(prop={'size':'small'},loc='best')
ax.grid(b=True,which='both',axis='x')
ax.grid(b=True,which='major',axis='y')

fig.tight_layout()
fig.savefig('vanGenuchten')

