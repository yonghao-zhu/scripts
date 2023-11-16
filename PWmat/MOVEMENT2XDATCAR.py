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
    convert from MOVEMENT (TDDFT pwmat) to XDATCAR (vasp)
    update: 2023.08.09
'''

import numpy as np
import linecache
import os

np.set_printoptions(precision=15)

#-----------------------------------------------------------------
def ReadMOVEMENT_NVE(fname_MOVEMENT, fname_POSCAR, run_type):
    # read natoms from POSCAR
    natoms_str = linecache.getline(fname_POSCAR, 7)
    natoms     = sum( [int(i) for i in natoms_str.split()] )
    print('--> natoms = ', natoms)

    if run_type == 'bomd':
        free_lines = 12; zeros_time_line = 14; num_ = 4
    elif run_type == 'tddft':
        free_lines = 5; zeros_time_line = 7; num_ = 3

    # read MOVEMENT
    if not os.path.isfile('movement.npy'):
        with open(fname_MOVEMENT, 'r') as f:
            num_lines = len( f.readlines() )
            times     = int( num_lines/(natoms*num_ + num_ + free_lines + 1) )
        print('--> times  = ', times)

        # save movement
        movement = np.zeros([times, natoms, 3], dtype=np.float64)
        for it in range(times):
            it_ini_lines = it*(natoms*num_ + num_ + free_lines + 1) + zeros_time_line
            if it % 1000 == 0:
                print('    --> it_ini_lines = ', it_ini_lines)
            for iatom in range(natoms):
                pos = linecache.getline(fname_MOVEMENT, it_ini_lines+iatom).split()
                movement[it, iatom, 0] = float(pos[1])
                movement[it, iatom, 1] = float(pos[2])
                movement[it, iatom, 2] = float(pos[3])
        
        np.save('movement.npy', movement)
    else:
        movement = np.load('movement.npy')

    return movement

#-----------------------------------------------------------------
def ReadMOVEMENT_NPT(fname_MOVEMENT):

    natoms = int(linecache.getline(fname_MOVEMENT, 1).split()[0])
    print('--> natoms = ', natoms)

    # read MOVEMENT   
    with open(fname_MOVEMENT, 'r') as f:
        num_lines = len( f.readlines() )
        times     = int( num_lines/(natoms*4 + 4 + 19 + 1) )
    print('--> times  = ', times)

    # save movement and lattice
    lattice  = np.zeros([times, 3, 3], dtype=np.float64)
    movement = np.zeros([times, natoms, 3], dtype=np.float64)

    for it in range(times):
        it_ini_lines = it*(natoms*4 + 4 + 19 + 1) + 12
        if it % 1000 == 0:
            print('    --> it_ini_lines = ', it_ini_lines)
        for iatom in range(natoms):
            pos = linecache.getline(fname_MOVEMENT, it_ini_lines+iatom+9).split()
            movement[it, iatom, 0] = float(pos[1])
            movement[it, iatom, 1] = float(pos[2])
            movement[it, iatom, 2] = float(pos[3])

            for i in range(1,4):
                lat = linecache.getline(fname_MOVEMENT, it_ini_lines+i).split()
                
                lattice[it, i-1, 0] = float(lat[0])
                lattice[it, i-1, 1] = float(lat[1])
                lattice[it, i-1, 2] = float(lat[2])

    np.save('movement.npy', movement)   
    np.save('lattice.npy', lattice)      

    return 1, 2
#-----------------------------------------------------------------
def SaveXDATCAR(fname_POSCAR, movement):

    # read natoms from POSCAR
    lattice_str = [linecache.getline(fname_POSCAR, i) for i in range(1, 8)]

    # x = 3
    times, natoms, x = np.shape(movement)

    # write xdatcar
    with open('XDATCAR', 'w+') as xdat:
        for i in lattice_str:
            xdat.writelines(i)

        for it in range(times):
            if it % 1000 == 0:
                print('    --> it = ', it)
                
            xdat.writelines('Direct configuration= '+str(it+1)+'\n')
            for iatm in range(natoms):
                pos = movement[it, iatm]
                xdat.writelines('     ' + str(pos[0]) + '     ')
                xdat.writelines(str(pos[1]) + '     ')
                xdat.writelines(str(pos[2]) + '     '+'\n')

#-----------------------------------------------------------------
def main():

    fname_MOVEMENT = 'MOVEMENT'
    run_type       = 'bomd' # bomd or tddft
    ensemble       = 'NPT' # NVE or NPT
    fname_POSCAR   = 'CONTCAR.vasp'

    if ensemble == 'NVE':
        # read movement
        movement = ReadMOVEMENT_NVE(fname_MOVEMENT, fname_POSCAR, run_type)
        # save xdatcar
        SaveXDATCAR(fname_POSCAR, movement)
    elif ensemble == 'NPT':
        # read movement
        ReadMOVEMENT_NPT(fname_MOVEMENT)

if __name__ == "__main__":
    main()