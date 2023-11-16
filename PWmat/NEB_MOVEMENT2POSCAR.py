#! /usr/bin/python3
# -*- conding=UTF-8 -*-
#  .--,       .--,
# ( (  \.---./  ) )
#  '.__/o   o\__.'
#     {=  ^  =}
#      >  -  <
#     /       \
#    //       \\
#   //|   .   |\\
#   "'\       /'"_.-~^`'-.
#      \  _  /--'         `
#    ___)( )(___
#   (((__) (__)))    
 
import numpy as np
import os, linecache
#-------------------------------------------------------------------------
def GetPos(movement_file, Iteration, imag, atom):

    if not os.path.isfile(movement_file):
        exit('No %s! Exitting...' %movement_file)

    with open(movement_file) as f:
        lines      = f.readlines()
        atoms_line = []
        for i in range(len(lines)):
            if 'atoms' in lines[i]:
                ite = int(lines[i].split()[3])
                if ite == Iteration:
                    atoms_line.append(i)
        atoms_line = atoms_line[-1*imag:]

    latt = np.zeros([imag, 3, 3], dtype=float)
    pos  = np.zeros([imag, sum(atom), 3], dtype=float)
    for i in range(imag):
        for j in range(3):
            latt_line   = linecache.getline(movement_file, atoms_line[i] + j + 3).split()
            latt[i,j,0] = float(latt_line[0])       
            latt[i,j,1] = float(latt_line[1])
            latt[i,j,2] = float(latt_line[2])       
        for p in range(sum(atom)):
            pos_line   = linecache.getline(movement_file, atoms_line[i] + p + 7).split()
            pos[i,p,0] = float(pos_line[1])
            pos[i,p,1] = float(pos_line[2])
            pos[i,p,2] = float(pos_line[3])

    return latt, pos

def WritePOSCAR(movement_file, Iteration, imag, element, atom):
    latt, pos = GetPos(movement_file, Iteration, imag, atom)    

    for i in range(imag):
        pos_file = 'POSCAR-%s.vasp' %i
        with open(pos_file, 'w+') as f:   
            f.writelines('NEB'+'\n'+'1.0'+'\n')
            for m in range(3):
                f.writelines(str(latt[i,m,0])+' '+str(latt[i,m,1])+' '+str(latt[i,m,2])+'\n')
            for j in element:
                f.writelines(j+'  ')
            f.writelines('\n')
            for k in atom:
                f.writelines(str(k)+'  ')
            f.writelines('\n'+'D'+'\n')                
            for p in range(sum(atom)):
                f.writelines(str(pos[i,p,0])+' '+str(pos[i,p,1])+' '+str(pos[i,p,2])+'\n')

#-------------------------------------------------------------------------
def main():
    
    movement_file = './MOVEMENT'
    Iteration     = 14 # final interation 
    imag          = 9  # total 
    element       = ['Hf', 'O']
    atom          = [4, 8]

    WritePOSCAR(movement_file, Iteration, imag, element, atom)

if __name__ == '__main__':
    main()
