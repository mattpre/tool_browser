# -*- coding: cp1252 -*-
import numpy as np
import math
import cPickle as pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import xml.etree.ElementTree as ET
from io import BytesIO

ET.register_namespace("", "http://www.w3.org/2000/svg")
        
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)

nx = 10
ny = 10
patches = []
cvect = []
coords = []
for kx in range(nx):
    for ky in range(ny):
        crds = np.array([[kx,kx,kx+1,kx+1],
                         [ky,ky+1,ky+1,ky]]).T
        cvect.append(kx+ky*2)
        poly = Polygon(crds)
        coords.append(crds)
        patches.append(poly)

pc = PatchCollection(patches,edgecolors=('none',))
pc.set_array(np.array(cvect))
ax.add_collection(pc)
cb = fig.colorbar(pc)
cb.set_label('colors',size=16)
cb.ax.tick_params(labelsize=16)

fig.canvas.draw()
colors = pc.get_facecolors()
pc.remove()
for kp in range(len(colors)):
    poly = Polygon(coords[kp],fc=colors[kp],ec='none')
    ax.add_patch(poly)
    an = ax.annotate('cos %i'%(cvect[kp]),xy=np.mean(coords[kp],0),ha='center',va='center')
    poly.set_gid('mypatch_{:05d}'.format(kp))
    an.set_gid('mytooltip_{:05d}'.format(kp))
        
ax.set_xlim(0,nx)
ax.set_ylim(0,ny)

fig.tight_layout()
f = BytesIO()
fig.savefig(f,format='svg')
plt.close(fig)

# --- Add interactivity ---

# Create XML tree from the SVG file.
tree, xmlid = ET.XMLID(f.getvalue())
tree.set('onload', 'init(evt)')

for i in patches:
    # Get the index of the shape
    index = patches.index(i)
    try:
        # Hide the tooltips
        tooltip = xmlid['mytooltip_{:05d}'.format(index)]
        tooltip.set('visibility', 'hidden')
        # Assign onmouseover and onmouseout callbacks to patches.
        mypatch = xmlid['mypatch_{:05d}'.format(index)]
        mypatch.set('onmouseover', "ShowTooltip(this)")
        mypatch.set('onmouseout', "HideTooltip(this)")
    except:
        pass

# This is the script defining the ShowTooltip and HideTooltip functions.
script = """
    <script type="text/ecmascript">
    <![CDATA[

    function init(evt) {
        if ( window.svgDocument == null ) {
            svgDocument = evt.target.ownerDocument;
            }
        }

    function ShowTooltip(obj) {
        var cur = obj.id.split("_")[1];
        var tip = svgDocument.getElementById('mytooltip_' + cur);
        tip.setAttribute('visibility',"visible");
        var patch = svgDocument.getElementById('mypatch_' + cur);
        patch.setAttribute('stroke',"red")
        }

    function HideTooltip(obj) {
        var cur = obj.id.split("_")[1];
        var tip = svgDocument.getElementById('mytooltip_' + cur);
        tip.setAttribute('visibility',"hidden");
        var patch = svgDocument.getElementById('mypatch_' + cur);
        patch.setAttribute('stroke',"none")
        }

    ]]>
    </script>
    """

# Insert the script at the top of the file and save it.
tree.insert(0, ET.XML(script))
ET.ElementTree(tree).write('fig.svg')

