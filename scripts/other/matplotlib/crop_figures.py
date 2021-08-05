from PIL import Image
import numpy as np



plist = []
for astr in ['A1','A2','B1','B2']:
    plist.append('MEC04_pl2526-9Gliss_supercrise_'+astr+'_v22')
    plist.append('MEC04_pl2526-9Gliss_2002_'+astr+'_v22')
    plist.append('MEC04_pl2526-9Gliss_Annee_moyenne_'+astr+'_v22')
plist.append('MEC04_pl2526-9Gliss_Annee_moyenne_v22_zonesv7_nu049partial')
plist.append('MEC04_pl2526-9Gliss_supercrise_v22a')
plist.append('MEC04_pl2526-9Gliss_2002_v22b')
plist.append('MEC04_pl2526-9Gliss_98-99_avecFlux_v23')
plist.append('MEC04_pl2526-9Gliss_98-99_v23')
plist.append('MEC04_pl2526-9Gliss_98-99_Wasserwert042_v23')
plist.append('MEC04_pl2526-9Gliss_98-99_A1_avecFlux_v23')
plist.append('MEC04_pl2526-9Gliss_98-99_A1_v23')
plist.append('MEC04_pl2526-9Gliss_98-99_B2_avecFlux_v23')
plist.append('MEC04_pl2526-9Gliss_98-99_B2_v23')

for f in plist:
    image=Image.open('disp_carte_t366_'+f+'.png')
    image.load()


    imageSize = image.size
    if imageSize[0]==1000 and imageSize[1]==1400:
        imageBox = image.getbbox()
        imageBox = (120,130,950,1280)
        cropped=image.crop(imageBox)
        cropped.save('disp_carte_t366_'+f+'.png')
        print f
