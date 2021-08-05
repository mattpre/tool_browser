# -*- coding: cp1252 -*-
import ZSmaterial,read_mat,re,sys
import xlsxwriter

def writeMat(inppath,outpath,columns):
    gd = ZSmaterial.GlobalData()

    materials = read_mat.read_mat(gd,inppath)

    workbook = xlsxwriter.Workbook(str(outpath[0]))

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
            for key in data.keys():
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
            try:
                unit = gd.unit_dict[head]
            except:
                unit = '[?]'
            try:
                astr = re.sub('L',gd.units[1],unit)
                astr = re.sub('F',gd.units[0],astr)
                astr = re.sub('D',u'°',astr)
                astr = re.sub('T',gd.units[3],astr)
                astr = re.sub('H',gd.units[4],astr)
                astr = re.sub('P',gd.units[5],astr)
            except:
                astr = unit
            worksheet.write_string(1,4+kh,astr,header)

        for km,mat in enumerate(materials):
##            name = unicode(mat.name,sys.stdin.encoding)
            name = mat.name
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
    if 'g' in columns:
        worksheet.write_rich_string(0,2+columns.index('g'),symbol,'g',header)
        worksheet.write_rich_string(1,2+columns.index('g'),'[kN/m',superscript,'3',']',header)
    if 'gD' in columns:
        worksheet.write_rich_string(0,2+columns.index('gD'),symbol,'g',subscript,'D',header)
        worksheet.write_rich_string(1,2+columns.index('gD'),'[kN/m',superscript,'3',']',header)
    if 'E50' in columns:
        worksheet.write_rich_string(0,2+columns.index('E50'),'E',subscript,'50',header)
        worksheet.write(1,2+columns.index('E50'),'[MPa]',header)
    if 'Eur' in columns:
        worksheet.write_rich_string(0,2+columns.index('Eur'),'E',subscript,'ur',header)
        worksheet.write(1,2+columns.index('Eur'),'[MPa]',header)
    if 'E0' in columns:
        worksheet.write_rich_string(0,2+columns.index('E0'),'E',subscript,'0',header)
        worksheet.write(1,2+columns.index('E0'),'[MPa]',header)
    if 'sref' in columns:
        worksheet.write_rich_string(0,2+columns.index('sref'),symbol,'s',subscript,'h,ref',header)
        worksheet.write(1,2+columns.index('sref'),'[kPa]',header)
    if 'phi' in columns:
        worksheet.write_rich_string(0,2+columns.index('phi'),symbol,'f','\'',header)
        worksheet.write(1,2+columns.index('phi'),u'[°]',header)
    if 'psi' in columns:
        worksheet.write_rich_string(0,2+columns.index('psi'),symbol,'y','\'',header)
        worksheet.write(1,2+columns.index('psi'),u'[°]',header)
    if 'c' in columns:
        worksheet.write(0,2+columns.index('c'),'c\'',header)
        worksheet.write(1,2+columns.index('c'),'[kPa]',header)
    if 'k_x' in columns:
        worksheet.write_rich_string(0,2+columns.index('k_x'),'k',subscript,'x',header)
        worksheet.write(1,2+columns.index('k_x'),'[m/s]',header)
    if 'k_y' in columns:
        worksheet.write_rich_string(0,2+columns.index('k_y'),'k',subscript,'y',header)
        worksheet.write(1,2+columns.index('k_y'),'[m/s]',header)
    if 'k_z' in columns:
        worksheet.write_rich_string(0,2+columns.index('k_z'),'k',subscript,'z',header)
        worksheet.write(1,2+columns.index('k_z'),'[m/s]',header)
    worksheet.set_column(0,0,12)
    worksheet.set_column(1,1,22)
    worksheet.set_column(2,12,8)

    # units:

    # write material data:
    kkm = -1
    for km,mat in enumerate(materials):
        if (mat.type=='HS-small strain stiffness'
            or mat.type=='Densification model'
            or 'Mohr' in mat.type
            or 'Hoek-Brown' in mat.type):
            print('%i'%(km),mat.type,mat.name)
            kkm += 1
##            worksheet.write(kkm+2,0, unicode(mat.name,sys.stdin.encoding),cell)
            worksheet.write(kkm+2,0, mat.name,cell)
            worksheet.write(kkm+2,1, mat.type,cell)
            if 'g' in columns:
                worksheet.write(kkm+2,2+columns.index('g'), mat.dens['gamma'],cell)
            if 'gD' in columns:
                try:
                    worksheet.write(kkm+2,2+columns.index('gD'), mat.dens['gamma_D'],cell)
                except:
                    pass
            if (mat.type=='HS-small strain stiffness' or
                mat.type=='Densification model'):
                if 'E50' in columns:
                    worksheet.write(kkm+2,2+columns.index('E50'), mat.nonl['Eref_50']*1e-3,cell)
                if 'Eur' in columns:
                    worksheet.write(kkm+2,2+columns.index('Eur'), mat.elas['Eref_ur']*1e-3,cell)
                if 'E0' in columns:
                    worksheet.write(kkm+2,2+columns.index('E0'), mat.elas['Eref_0']*1e-3,cell)
                if 'sref' in columns:
                    worksheet.write(kkm+2,2+columns.index('sref'), mat.elas['sig_ref'],cell)
            else:
                if 'E50' in columns and 'E0' in columns:
                    worksheet.merge_range(kkm+2,2+columns.index('E50'),kkm+2,2+columns.index('E0'), mat.elas['E']*1e-3,cell)
                elif 'E50' in columns and 'Eur' in columns:
                    worksheet.merge_range(kkm+2,2+columns.index('E50'),kkm+2,2+columns.index('Eur'), mat.elas['E']*1e-3,cell)
            try:
                if 'phi' in columns:
                    worksheet.write(kkm+2,2+columns.index('phi'), mat.nonl['phi'],cell)
                if 'psi' in columns:
                    worksheet.write(kkm+2,2+columns.index('psi'), mat.nonl['psi'],cell)
                if 'c' in columns:
                    worksheet.write(kkm+2,2+columns.index('c'), mat.nonl['c'],cell)
            except:
                if 'phi' in columns:
                    worksheet.write(kkm+2,2+columns.index('phi'),'',cell)
                if 'psi' in columns:
                    worksheet.write(kkm+2,2+columns.index('psi'),'',cell)
                if 'c' in columns:
                    worksheet.write(kkm+2,2+columns.index('c'),'',cell)
            try:
                if 'k_x' in columns:
                    worksheet.write(kkm+2,2+columns.index('k_x'), mat.flow['k_x'],cell)
            except:
                if 'k_x' in columns:
                    worksheet.write(kkm+2,2+columns.index('k_x'),'',cell)
            try:
                if 'k_y' in columns:
                    worksheet.write(kkm+2,2+columns.index('k_y'), mat.flow['k_y'],cell)
            except:
                if 'k_y' in columns:
                    worksheet.write(kkm+2,2+columns.index('k_y'),'',cell)
            try:
                if 'k_z' in columns:
                    worksheet.write(kkm+2,2+columns.index('k_z'), mat.flow['k_z'],cell)
            except:
                if 'k_z' in columns:
                    worksheet.write(kkm+2,2+columns.index('k_z'),'',cell)

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

