# -*- coding: cp1252 -*-

class GlobalData:
    def __init__(self):
        self.units = []
        self.nMat = 0
        self.unit_dict = {'M':'[-]',
                          'K0y':'[-]',
                          'K0x':'[-]',
                          'gamma_D':'[F/L3]',
                          'KSR_0':'[-]',
                          'nu_ur':'[-]',
                          'OCR':'[-]',
                          'Eref_ur':'[P]',
                          'Eref_50':'[P]',
                          'e_0':'[-]',
                          'gamma_f':'[F/L3]',
                          'INT5':'[-]',
                          'INT4':'[-]',
                          'INT7':'[-]',
                          'INT6':'[-]',
                          'INT1':'[-]',
                          'INT3':'[-]',
                          'E_oed':'[P]',
                          'Eref_0':'[P]',
                          'INT9':'[-]',
                          'nu':'[-]',
                          'phi':'[D]',
                          'psi':'[D]',
                          'E':'[P]',
                          'D':'[-]',
                          'H':'[-]',
                          'INT2':'[-]',
                          'SS_ext':'[-]',
                          'pmin_co':'[P]',
                          'sigref_oed':'[P]',
                          'R_f':'[-]',
                          'c':'[P]',
                          'sig_L':'[P]',
                          'KNC_0':'[-]',
                          'm':'[-]',
                          'gam_07':'[-]',
                          'INT8':'[-]',
                          'f_t':'[P]',
                          'f_c':'[P]',
                          'sig_ref':'[P]',
                          'gamma':'[F/L3]',
                          'Ix':'[L4]',
                          'Iy':'[L4]',
                          'Iz':'[L4]',
                          'b':'[L]',
                          'h':'[L]',
                          'inp_type':'[-]',
                          'Ax':'[L2]',
                          'Ay':'[L2]',
                          'Az':'[L2]',
                          'Kt/Kn':'[-]',
                          'Kn':'[-]',
                          'inherit':'[-]',
                          'k_x':'[L/T]',
                          'k_y':'[L/T]',
                          'k_z':'[L/T]',
                          'S_rres':'[-]',
                          'alpha':'[-]',
                          'pen':'[-]',
                          'cutoff':'[-]',
                          'K_f':'[-]',
                          'Profile':'[-]',
                          'Database':'[-]',
                          'm_b':'[-]',
                          'a':'[-]',
                          's':'[-]'}

class ZSmaterial:
    def __init__(self):
        self.version = 0
        self.number = 0
        self.name = 0
        self.type = 0
        self.buttons = []
        self.groups = []
        
        self.elas = {}
        self.main = {}
        self.nonl = {}
        self.geom = {}
        self.dens = {}
        self.flow = {}
        self.creep = {}
        self.heat = {}
        self.humid = {}
        self.inis = {}
        self.stab = {}
        self.damp = {}

    def read_Data(self,data):
        continuum = ['Elastic','HS-small strain stiffness','Mohr-Coulomb']
        
        lines = data[0]
        self.groups.append('elas')
        if 'Elastic' in self.type or 'Mohr-Coulomb' in self.type \
           or self.type=='Multilaminate' or 'Hoek-Brown' in self.type:
            self.elas['E'] = float(lines[0].split()[1])
            self.elas['nu'] = float(lines[0].split()[2])
        elif (self.type=='HS-small strain stiffness' or
              self.type=='Densification model'):
            self.elas['Eref_ur'] = float(lines[0].split()[1])
            self.elas['sig_ref'] = float(lines[0].split()[4])
            self.elas['nu_ur'] = float(lines[0].split()[2])
            self.elas['m'] = float(lines[0].split()[3])
            self.elas['sig_L'] = float(lines[0].split()[5])
            self.elas['SS_ext'] = int(float(lines[0].split()[7]))
            self.elas['Eref_0'] = float(lines[0].split()[8])
            self.elas['gam_07'] = float(lines[0].split()[9])
        elif self.type=='Contact':
            self.elas['Kn'] = float(lines[0].split()[1])
            self.elas['Kt/Kn'] = float(lines[0].split()[2])

        lines = data[1]
        if self.type=='Elastic Beam':
            self.groups.append('geom')
            self.geom['inp_type'] = int(lines[0].split()[1])
            if self.geom['inp_type']==0:  # specification by profile DB
                self.geom['Database'] = lines[1]
                self.geom['Profile'] = lines[2]
                ind = 3
            elif self.geom['inp_type']==1:  # specification by geometry
                self.geom['b'] = float(lines[1].split()[1])
                self.geom['h'] = float(lines[1].split()[3])
                ind = 2
            elif self.geom['inp_type']==2:  # specification by values
                ind = 1
            else:
                print('Error in read_Data')
            self.geom['Ix'] = float(lines[ind].split()[3])
            self.geom['Iy'] = float(lines[ind].split()[4])
            self.geom['Iz'] = float(lines[ind].split()[5])
            self.geom['Ax'] = float(lines[ind].split()[0])
            self.geom['Ay'] = float(lines[ind].split()[1])
            self.geom['Az'] = float(lines[ind].split()[2])

        lines = data[2]
        if self.type=='Beam':
            self.groups.append('main')
            self.main['flex_based'] = int(lines[0].split()[1])

        lines = data[3]
        if 'interface' not in self.type and 'hinge' not in self.type:
            self.groups.append('dens')
            self.dens['gamma'] = float(lines[0].split()[1])
            if self.type in continuum:
                self.dens['gamma_f'] = float(lines[0].split()[2])
                self.dens['e_0'] = float(lines[0].split()[3])
                self.dens['gamma_D'] = float(lines[0].split()[6])

        lines = data[4]
        if self.buttons[4]==1 and not self.type=='Seepage' and 'hinge' not in self.type:
            self.groups.append('flow')
            self.flow['k_x'] = float(lines[0].split()[1])
            self.flow['k_y'] = float(lines[0].split()[2])
            self.flow['k_z'] = float(lines[0].split()[3])
            self.flow['K_f'] = float(lines[0].split()[11])
            if self.type in continuum:
                self.flow['S_rres'] = float(lines[0].split()[12])
                self.flow['alpha'] = float(lines[0].split()[13])
                self.flow['pen'] = float(lines[0].split()[17])
                self.flow['cutoff'] = float(lines[0].split()[18])

        lines = data[5]
        if self.buttons[5]==1:
            self.groups.append('creep')
            self.creep['Ad'] = float(lines[0].split()[4])
            self.creep['Av'] = float(lines[0].split()[6])
            self.creep['Bd'] = float(lines[0].split()[5])
            self.creep['Bv'] = float(lines[0].split()[7])
            self.creep['EF1'] = int(lines[0].split()[2])
            self.creep['EF2'] = int(lines[0].split()[3])
            self.creep['a'] = float(lines[0].split()[8])
            self.creep['b'] = float(lines[0].split()[9])

        lines = data[6]
        if 'Mohr-Coulomb' in self.type:
            self.groups.append('nonl')
            self.nonl['c'] = float(lines[0].split()[4])
            self.nonl['phi'] = float(lines[0].split()[2])
            self.nonl['psi'] = float(lines[0].split()[3])
        elif 'Hoek-Brown' in self.type:
            self.groups.append('nonl')
            self.nonl['f_c'] = float(lines[1].split()[0])
##            self.nonl['psi'] = float(lines[0].split()[3])
            self.nonl['m'] = float(lines[1].split()[1])
            self.nonl['s'] = float(lines[1].split()[2])
            self.nonl['a'] = float(lines[1].split()[3])
        elif (self.type=='HS-small strain stiffness' or
              self.type=='Densification model'):
            self.groups.append('nonl')
            self.nonl['phi'] = float(lines[0].split()[1])
            self.nonl['psi'] = float(lines[0].split()[2])
            self.nonl['c'] = float(lines[0].split()[3])
            self.nonl['Eref_50'] = float(lines[0].split()[4])
            self.nonl['R_f'] = float(lines[0].split()[5])
            self.nonl['f_t'] = int(float(lines[0].split()[6]))
            self.nonl['D'] = float(lines[0].split()[7])
            self.nonl['INT1'] = int(float(lines[0].split()[8]))
            self.nonl['INT2'] = int(float(lines[0].split()[9]))
            self.nonl['INT3'] = int(float(lines[0].split()[10]))
            self.nonl['H'] = float(lines[1].split()[0])
            self.nonl['M'] = float(lines[1].split()[1])
            self.nonl['KNC_0'] = float(lines[1].split()[2])
            self.nonl['sigref_oed'] = float(lines[1].split()[3])
            self.nonl['E_oed'] = float(lines[1].split()[4])
            self.nonl['INT4'] = int(float(lines[2].split()[0]))
            self.nonl['INT5'] = int(float(lines[2].split()[1]))
            self.nonl['INT6'] = int(float(lines[3].split()[0]))
            self.nonl['pmin_co'] = int(float(lines[3].split()[1]))
            self.nonl['INT7'] = int(float(lines[3].split()[2]))
            self.nonl['OCR'] = int(float(lines[3].split()[3]))
            self.nonl['KSR_0'] = int(float(lines[3].split()[4]))
            self.nonl['INT8'] = int(float(lines[3].split()[5]))
            self.nonl['INT9'] = int(float(lines[3].split()[6]))
        elif self.type=='Contact':
            self.groups.append('nonl')
            if False:   # pre-v20 format?
                self.nonl['inherit'] = float(lines[0].split()[5])
                if self.nonl['inherit']==0:
                    self.nonl['phi'] = float(lines[0].split()[1])
                    self.nonl['psi'] = float(lines[0].split()[2])
                    self.nonl['c'] = float(lines[0].split()[3])
            else:
                self.nonl['inherit'] = float(lines[1].split()[3])
                if self.nonl['inherit']==0:
                    self.nonl['phi'] = float(lines[1].split()[0])
                    self.nonl['psi'] = float(lines[1].split()[1])
                    self.nonl['c'] = float(lines[1].split()[2])

        lines = data[7]
        if self.buttons[7]==1:
            self.groups.append('heat')
            self.heat['dil'] = float(lines[0].split()[1])
            if self.type in continuum:
                self.heat['cond'] = float(lines[0].split()[2])
                self.heat['cap'] = float(lines[0].split()[3])

        lines = data[8]
        if self.buttons[8]==1:
            self.groups.append('humid')
            self.humid['dil'] = float(lines[0].split()[1])
            self.humid['a'] = float(lines[0].split()[2])
            self.humid['WI'] = float(lines[0].split()[3])
            self.humid['DI'] = float(lines[0].split()[4])

        lines = data[9]
        if self.buttons[9]==1:
            self.groups.append('inis')
            self.inis['K0x'] = float(lines[0].split()[1])
            self.inis['K0y'] = float(lines[0].split()[2])

        lines = data[12]
        if self.buttons[13]==1:
            self.groups.append('damp')
            self.damp['alpha'] = float(lines[0].split()[1])
            self.damp['beta'] = float(lines[0].split()[2])
    
