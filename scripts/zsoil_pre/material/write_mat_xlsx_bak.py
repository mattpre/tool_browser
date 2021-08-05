# -*- coding: cp1252 -*-
import ZSmaterial,read_mat,re,sys
import xlsxwriter

gd = ZSmaterial.GlobalData()

path = '//CALCUL-FEV13/Mandats sur F RAID0/M843_CEVA_sortie3'
prob = 'M843_sortie3_etapes_piedroitNoCnt'
path = '.'
prob = 'test_tube'
path = 'D:/Mandats/M921_MurTuescherz/technique/fev2015'
prob = 'M921_fev2015_v1_d15e4'
##path = '//CALCUL-APR11/Mandats Disc C/M877_Braunwald'
##prob = 'MEC04_pl2526-9Gliss_2002_A2_v22'
##path = '//192.168.1.53/Mandats_53_E/M883_Zug'
##prob = 'M883_Alpenstrasse_Grundbruch'
##path = '//CALCUL-FEV13/Mandats sur F RAID0/M972_Frank_Thomas'
##prob = 'M972_FT_mat_v3'
path = 'D:/Mandats/M971_Vauseyon/technique/2D'
prob = 'M971_VM1_v4_macNL_marneML'
path = '//CALCUL-APR11/Mandats Disc C/M982_LEB'
prob = 'M982_PROFIL4_3D_V3_GEOL1_CNT'
path = 'D:/Mandats/M985_Chaudanne/technique'
prob = 'M985'
path = 'D:/Mandats/M630_CEVA/SOVACB/technique/calage'
prob = 'M630_profil47_EurE0'
path = 'D:/Mandats/M992_Ottikerstrasse/technique'
prob = 'M992_2D'
path = 'D:/Mandats/M1004_Menuisiers/technique'
prob = 'M1004_pieu_axi_geoBase'
##path = '//CALCUL-APR11/Mandats Disc D/M911_SOVALP_Lancy/structure'
##prob = 'M911_Sovalp_3D-v1_7_structure_LpieuxParams'
##path = 'D:/Mandats/M972_CEVA_FrankThomas/technique/expertise_nov2015/profil10'
##prob = 'p10_nappe3b_f388_5_Eb30_sradl50_6c2-6c12_EurE0'
materials = read_mat.read_mat(gd,path + '/' + prob + '.inp')

workbook = xlsxwriter.Workbook(path + '/mat_' + prob + '.xlsx')

# format definitions:
subscript = workbook.add_format({'font_script':2})
superscript = workbook.add_format({'font_script':1})
symbol = workbook.add_format()
symbol.set_font_name('Symbol')
header = workbook.add_format({'bg_color':'#C5BE97','border':1,'align':'vcenter'})

groupnames = ['Elastic','Geometry','Main','Density','Flow',
              'Creep','Nonlinear','Heat','Humidity','Stability','Damping',
              'Initial state']
groupnames_active = []
gps = ['elas','geom','main','dens','flow',
       'creep','nonl','heat','humid','stab','damp','inis']

gps_active0 = set()
gps_active = []
for km,mat in enumerate(materials):
    for gp in mat.groups:
        gps_active0.add(gp)
for kgp,gp in enumerate(gps):
    if gp in gps_active0:
        gps_active.append(gp)
        groupnames_active.append(groupnames[kgp])
nGroups = len(gps_active)

headers = [set() for kg in range(nGroups)]
for km,mat in enumerate(materials):
    for kg in range(nGroups):
        data = getattr(mat,gps_active[kg])
        for key in data.iterkeys():
            headers[kg].add(key)
for kg,h in enumerate(headers):
    headers[kg] = list(h)

for kg in range(nGroups):
    worksheet = workbook.add_worksheet(groupnames_active[kg])

    # write general headers:
    worksheet.write(0,0,'Version',header)
    worksheet.write(0,1,'Number',header)
    worksheet.write(0,2,'Name',header)
    worksheet.write(0,3,'Type',header)
    worksheet.set_column(0,0,6.86)
    worksheet.set_column(1,1,7.29)
    worksheet.set_column(2,2,30)
    worksheet.set_column(3,3,21.86)

    # write group-specific headers:
    for kh,head in enumerate(headers[kg]):
        if '_' in head:
            worksheet.write_rich_string(0,4+kh,head.split('_')[0],subscript,head.split('_')[1],header)
        else:
            worksheet.write_rich_string(0,4+kh,head,header)

    # write units:
    for kh,head in enumerate(headers[kg]):
        astr = re.sub('L',gd.units[1],gd.unit_dict[head])
        astr = re.sub('F',gd.units[0],astr)
        astr = re.sub('D',u'°',astr)
        astr = re.sub('T',gd.units[3],astr)
        astr = re.sub('H',gd.units[4],astr)
        astr = re.sub('P',gd.units[5],astr)
        worksheet.write_string(1,4+kh,astr,header)

    for km,mat in enumerate(materials):
        name = unicode(mat.name,sys.stdin.encoding)
        worksheet.write(km+2,0, 'v'+mat.version)
        worksheet.write_number(km+2,1, int(mat.number))
        worksheet.write(km+2,2, name)
        worksheet.write(km+2,3, mat.type)

        data = getattr(mat,gps_active[kg])
        for key in data:
            worksheet.write(km+2,4+headers[kg].index(key),data[key])

# create summary continuum:

header = workbook.add_format({'bg_color':'#C5BE97','border':1,'align':'center'})
header.set_align('vcenter')
cell = workbook.add_format({'border':1,'align':'center'})
cell.set_align('vcenter')

worksheet = workbook.add_worksheet('Summary')

# headers:
worksheet.merge_range(0,0,1,0,u'Matériau',header)
worksheet.merge_range(0,1,1,1,'Type',header)
worksheet.write_rich_string(0,2,symbol,'g',header)
worksheet.write_rich_string(0,3,symbol,'g',subscript,'D',header)
worksheet.write_rich_string(0,4,'E',subscript,'50',header)
worksheet.write_rich_string(0,5,'E',subscript,'ur',header)
worksheet.write_rich_string(0,6,'E',subscript,'0',header)
worksheet.write_rich_string(0,7,symbol,'s',subscript,'h,ref',header)
worksheet.write_rich_string(0,8,symbol,'f','\'',header)
worksheet.write(0,9,'c\'',header)
worksheet.set_column(0,0,12)
worksheet.set_column(1,1,22)
worksheet.set_column(2,12,8)
##worksheet.set_column(3,3,5)
##worksheet.set_column(3,3,5)
##worksheet.set_column(3,3,5)
##worksheet.set_column(3,3,5)
##worksheet.set_column(3,3,5)

# units:
worksheet.write_rich_string(1,2,'[kN/m',superscript,'3',']',header)
worksheet.write_rich_string(1,3,'[kN/m',superscript,'3',']',header)
worksheet.write(1,4,'[MPa]',header)
worksheet.write(1,5,'[MPa]',header)
worksheet.write(1,6,'[MPa]',header)
worksheet.write(1,7,'[kPa]',header)
worksheet.write(1,8,u'[°]',header)
worksheet.write(1,9,'[kPa]',header)

# write material data:
kkm = -1
for km,mat in enumerate(materials):
    if (mat.type=='HS-small strain stiffness'
        or mat.type=='Densification model'
        or 'Mohr' in mat.type
        or 'Hoek-Brown' in mat.type):
        print km,mat.type,mat.name
        kkm += 1
        worksheet.write(kkm+2,0, unicode(mat.name,sys.stdin.encoding),cell)
        worksheet.write(kkm+2,1, mat.type,cell)
        worksheet.write(kkm+2,2, mat.dens['gamma'],cell)
        try:
            worksheet.write(kkm+2,3, mat.dens['gamma_D'],cell)
        except:
            pass
        if (mat.type=='HS-small strain stiffness' or
            mat.type=='Densification model'):
            worksheet.write(kkm+2,4, mat.nonl['Eref_50']*1e-3,cell)
            worksheet.write(kkm+2,5, mat.elas['Eref_ur']*1e-3,cell)
            worksheet.write(kkm+2,6, mat.elas['Eref_0']*1e-3,cell)
            worksheet.write(kkm+2,7, mat.elas['sig_ref'],cell)
        else:
            worksheet.merge_range(kkm+2,4,kkm+2,6, mat.elas['E']*1e-3,cell)
        try:
            worksheet.write(kkm+2,8, mat.nonl['phi'],cell)
            worksheet.write(kkm+2,9, mat.nonl['c'],cell)
        except:
            worksheet.write(kkm+2,8, '',cell)
            worksheet.write(kkm+2,9, '',cell)
        

workbook.close()

### create summary reste:
##
##header = workbook.add_format({'bg_color':'#C5BE97','border':1,'align':'center'})
##header.set_align('vcenter')
##cell = workbook.add_format({'border':1,'align':'center'})
##cell.set_align('vcenter')
##
##worksheet = workbook.add_worksheet('Summary')
##
### headers:
##worksheet.merge_range(0,0,1,0,u'Matériau',header)
##worksheet.merge_range(0,1,1,1,'Type',header)
##worksheet.write(0,2,u'Densité',header)
##worksheet.write_rich_string(0,3,'E',subscript,'50',header)
##worksheet.write_rich_string(0,4,'E',subscript,'ur',header)
##worksheet.write_rich_string(0,5,'E',subscript,'0',header)
##worksheet.write_rich_string(0,6,symbol,'f','\'',header)
##worksheet.write(0,7,'c\'',header)
##worksheet.set_column(0,0,12)
##worksheet.set_column(1,1,15)
##worksheet.set_column(2,12,8)
##worksheet.set_column(2,12,6)
####worksheet.set_column(3,3,5)
####worksheet.set_column(3,3,5)
####worksheet.set_column(3,3,5)
####worksheet.set_column(3,3,5)
####worksheet.set_column(3,3,5)
##
### units:
##worksheet.write_rich_string(1,2,'[kN/m',superscript,'3',']',header)
##worksheet.write(1,3,'[MPa]',header)
##worksheet.write(1,4,'[MPa]',header)
##worksheet.write(1,5,'[MPa]',header)
##worksheet.write(1,6,u'[°]',header)
##worksheet.write(1,7,'[kPa]',header)
##
### write material data:
##kkm = -1
##for km,mat in enumerate(materials):
##    if not (mat.type=='Contact'
##            or mat.type=='Nail interface'
##            or 'Beam' in mat.type
##            or 'Truss' in mat.type
##            or 'Orthotropic' in mat.type
##            or 'Seepage' in mat.type):
##        print km,mat.type,mat.name
##        kkm += 1
##        worksheet.write(kkm+2,0, unicode(mat.name,sys.stdin.encoding),cell)
##        worksheet.write(kkm+2,1, mat.type,cell)
##        worksheet.write(kkm+2,2, mat.dens['gamma'],cell)
##        if (mat.type=='HS-small strain stiffness' or
##            mat.type=='Densification model'):
##            worksheet.write(kkm+2,3, mat.nonl['Eref_50']*1e-3,cell)
##            worksheet.write(kkm+2,4, mat.elas['Eref_ur']*1e-3,cell)
##            worksheet.write(kkm+2,5, mat.elas['Eref_0']*1e-3,cell)
##        else:
##            worksheet.merge_range(kkm+2,3,kkm+2,5, mat.elas['E']*1e-3,cell)
##        try:
##            worksheet.write(kkm+2,6, mat.nonl['phi'],cell)
##            worksheet.write(kkm+2,7, mat.nonl['c'],cell)
##        except:
##            worksheet.write(kkm+2,6, '',cell)
##            worksheet.write(kkm+2,7, '',cell)
##        
##
##workbook.close()

