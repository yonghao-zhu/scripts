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
import os
import linecache

import matplotlib as mpl
import matplotlib.pyplot as plt

#--------------------------------------------------------------------
def GetCenter(ensemble, movement, lattice, atoms_rang, time_range, direction, save_=True):

    if direction == 'x':
        direc = 0
    elif direction == 'y':
        direc = 1
    elif direction == 'z':
        direc = 2
    else:
        exit('Error direction!')

    atoms_rang = np.array(atoms_rang)
    time_range = np.array([i for i in range(time_range[0]-1, time_range[1], 1)])

    movement   = np.load(movement)[:, atoms_rang][time_range]
    print('movement.shape = ', movement.shape)    

    if ensemble == 'NPT':
        lattice = np.load(lattice)

    if len(np.where(movement[:] < 0.0)[0]) > 0 or len(np.where(movement[:] > 1.0)[0]) > 0:
        print('Warnning: exceed boundary!')

    if ensemble == 'NPT':
        for it in range(np.shape(time_range)[0]):
            movement[it] = np.dot(movement[it], lattice[it])

    center = np.mean(movement[:, direc], axis=1)

    if save_:
        np.savetxt('center-1.txt', center)

    return center

def PlotTime(center, time_range):

    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.18)

    plt.plot(time_range, center, lw=2)  

    plt.show()

#--------------------------------------------------------------------
def main():

    atoms_rang = [i for i in range(54, 108)] # list, start 0
    ensemble   = 'NPT' # NPT or NVE
    movement   = 'movement.npy' 
    lattice    = 'lattice.npy' 
    time_range = [1, 1000] # MD steps
    step       = 1 # fs
    direction  = 'z' # x, y, z

    # get time fs
    time_range_fs = [i*step for i in range(time_range[0], time_range[1]+1, 1)]
    print('total times    = ', len(time_range_fs), 'fs')

    center = GetCenter(ensemble, movement, lattice, atoms_rang, time_range, direction)

    PlotTime(center, time_range_fs)

if __name__ == '__main__':
    main()
