

from zsoil_tools import zsoil_inp as zi
import sys,math,numpy
import vtk



pathname = '//CALCUL-APR11/Mandats Disc C/M974_Chatelard'
pathname = '..'
prob = 'M1224_3Dcomplet_v1_7'

mesh = zi(pathname,prob)
mesh.read_inp()
mesh.write_vtu(prob,steps=[0])



        
