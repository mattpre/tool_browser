

import ZSmaterial


def read_mat(global_data,fname):    
    f = open(fname)

    materials = []
    for line in iter(lambda: f.readline(), ""):
        if 'NUM_MATERIALS=' in line:
            global_data.nMat = int(float(line.split('=')[1]))
            line = f.readline()
            for km in range(global_data.nMat):
                mat = ZSmaterial.ZSmaterial()
                materials.append(mat)
                # readmat proc:
                mat.version = line.split()[1]
                mat.number = line.split()[2]
                line = f.readline()
                mat.type = line[:-1]
                line = f.readline()
                mat.name = line[:-1]
                line = f.readline()
                v = line.split('=')[1].split()
                for val in v:
                    mat.buttons.append(int(float(val)))
                data = [[] for k in range(13)]
                line = f.readline()
                while not 'MATERIAL' in line:
                    if 'ELAS' in line:
                        ind = 0
                        data[0].append(line)
                    elif 'GEOM' in line:
                        ind = 1
                        data[1].append(line)
                    elif 'MAIN' in line:
                        ind = 2
                        data[2].append(line)
                    elif 'DENS' in line:
                        ind = 3
                        data[3].append(line)
                    elif 'FLOW' in line:
                        ind = 4
                        data[4].append(line)
                    elif 'CREEP' in line:
                        ind = 5
                        data[5].append(line)
                    elif 'NONL' in line:
                        ind = 6
                        data[6].append(line)
                    elif 'HEAT' in line:
                        ind = 7
                        data[7].append(line)
                    elif 'HUMID' in line:
                        ind = 8
                        data[8].append(line)
                    elif 'INIS' in line:
                        ind = 9
                        data[9].append(line)
                    elif 'STAB' in line:
                        ind = 10
                        data[10].append(line)
                    elif 'DISC' in line:
                        ind = 11
                        data[11].append(line)
                    elif 'DAMP' in line:
                        ind = 12
                        data[12].append(line)
                    elif not line=='\n':
                        data[ind].append(line)
                    line = f.readline()
                print(mat.number,mat.name,mat.type)
                mat.read_Data(data)
                # end readmat proc
        elif 'STANDARD' in line:
            line = f.readline()
            v = line.split()
            for vv in v:
                global_data.units.append(vv)
            global_data.units.append(v[0][0]+'Pa')

    f.close()

    return materials
