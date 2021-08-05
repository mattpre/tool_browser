

import os,sys
import shutil
import re

def pathsep(path):
    folders = []
    while 1:
        path, folder = os.path.split(path)

        if folder != "":
            folders.append(folder)
        else:
            if path != "":
                folders.append(path)
            break
    folders.reverse()
    return folders

def subdir(path,level):
    n = len(pathsep(path))
    if level>=n:
        return False
    else:
        subdir = ''
        for kk in range(n-level):
            path, folder = os.path.split(path)

            if folder != "":
                subdir = folder+'/'+subdir
            else:
                break
        return subdir[:-1]

extlist = ('.inp','.docx','.xlsx','.pdf','.pptx','.ppt','.xls','.doc','.py','.log','.png','.csv')

# get list of mandats:
Mandats = []
mdict = {}
for root in os.listdir('Z:'):
    if os.path.isdir('Z:/'+root):
        mandat = root.split('_')[0]
        if 'M' in mandat:
            if mandat[1:].isdigit():
                mNo = int(mandat[1:])
                if not mNo in Mandats:
                    Mandats.append(mNo)
                    mdict[mNo] = root

for root, dirs, files in os.walk(r'.'):
    if 'M1036' in root:
        print(root)
    if len(pathsep(root))>3:
        for file in files:
            path_file = os.path.join(root,file)
            parts = pathsep(root)
            pind = 0
            for kp in range(len(parts)):
                if re.match('M[19][0-9]{2,3}',parts[kp]):
                    pind = kp
                    break
            if file.endswith(extlist):
                mandat = parts[pind].split('_')[0]
                if 'M' in mandat:
                    if mandat[1:].isdigit():
                        mNo = int(mandat[1:])
                        if not mNo in Mandats:
                            print('%s does not exist'%(parts[pind]))
                            newroot = parts[pind]
                            os.mkdir('Z:/'+newroot)
                            Mandats.append(mNo)
                            mdict[mNo] = newroot
                        else:
                            newroot = mdict[mNo]+'/51'
                        if len(parts)>pind+1:
                            newpath = 'Z:/'+newroot+'/'+subdir(root,pind+1)
                        else:
                            newpath = 'Z:/'+newroot
##                            print(newpath+'/'+file)
##                        open(path_file).close()
##                        open('Z:/'+root.split('\\',5)[-1]+'/'+file).close()
                        if not os.path.isdir(newpath):
                            print(newpath+'  path created')
                            os.makedirs(newpath)
                        if not os.path.isfile(newpath+'/'+file):
                            if os.path.isfile(path_file):
                                print('new: '+newpath+'/'+file)
                                shutil.copy2(path_file,newpath+'/'+file)
                        else:
                            if os.path.getmtime(path_file)>os.path.getmtime(newpath+'/'+file):
                                print('rep: '+newpath+'/'+file)
                                shutil.copy2(path_file,newpath+'/'+file)
                            else:
                                print('kep: '+newpath+'/'+file)
##                        quit()
##    except:
##        print(file)
##        pass
##                    print(int(mandat[1:]),file)
##        shutil.copy2(path_file,'destination_directory') # change you destination dir


