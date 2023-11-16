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
 
'''
    plot the TDDOS
    need RunDos/dos-xxxx/charge/RHO.x-x.xsf dir
'''

import numpy as np
np.set_printoptions(precision=15)
import linecache
import os

#=========================================================================
def MoveAtomsFromPoscar(poscar, move_dis, move_atom, save_=True):

    if not os.path.isfile(poscar):
        exit('No POSCAR! Exiting...')
    if len(move_dis) != len(move_atom):
        exit('Error Inputs (len(move_dis) != len(move_atom))! Exitiing')

    # read poscar
    head = [linecache.getline(poscar, i) for i in range(8)]
    mode = linecache.getline(poscar, 8)
    atom = sum([int(i) for i in linecache.getline(poscar, 7).split()])
    latt = np.zeros([3,3], dtype=float)
    for i in range(3):
        tmp = linecache.getline(poscar, 3+i).split()
        for j in range(3):
            latt[i,j] = float(tmp[j])
    print('Total Atoms = ', atom) 
    if 'D' not in mode:
        exit('Please input Direct mode of POSCAR!')
    pos = np.zeros([atom, 3], dtype=float)
    for ia in range(atom):
        tmp = linecache.getline(poscar, 9 + ia).split()
        pos[ia, 0] = float(tmp[0])
        pos[ia, 1] = float(tmp[1])
        pos[ia, 2] = float(tmp[2])

    # move atom
    pos = np.dot(pos, latt)
    for i in range(len(move_dis)):
        iatom = np.array(move_atom[i])
        pos[iatom] += np.array(move_dis[i])

    # save
    if save_:
        with open('POSCAR.new.vasp', 'w+') as f:
            for i in head:
                f.writelines(i)
            f.writelines('C'+'\n')
            for ia in range(atom):
                for i in range(3):
                    f.writelines('    '+str(pos[ia,i]))
                f.writelines('\n')
                
#=========================================================================
def main():

    poscar    = 'POSCAR.vasp'
    move_dis  = [[0.0, 0.0, -0.3292761611363999/2], [0.0, 0.0, 0.3292761611363999/2]]
    # Ga and N
    move_atom = [[i for i in range(96)], [i for i in range(96,192)]]
    MoveAtomsFromPoscar(poscar, move_dis, move_atom)

if __name__ == '__main__':
    main()

