# @description Moves nodal coordinates
# @input inp-file
# @output inp-file with new '.ing'-section
# @mandat M1013 Hardturm
# @author Matthias Preisig


f = open('M1013_W_v6_pilesV1bis2.inp')
of = open('M1013_W_v6_pilesV1bis2-2.inp','w')

for line in f:
    if '.ing' in line:
        of.write(line)
        line = f.next()
        while len(line)>2:
            v = line.split()
            x = float(v[1])
            z = float(v[3])
            of.write('%s %1.8g %s %1.8g 0\n'%\
                     (v[0],x+257.37+130.85,v[2],z++323.65-317.15))
            line = f.next()
    else:
        of.write(line)    

of.close()
f.close()
