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
    generate the IN.TDDFT_TIME for TDDFT calculations
'''

import numpy as np
np.set_printoptions(precision=15)

import matplotlib as mpl
import matplotlib.pyplot as plt

#--------------------------------------------------------------------
'''
    ftddft = strength * time
'''
def Gen_LinearField(time, step, strength, save_=True):

    IN_TDDFT_TIME = np.zeros([int(time/step)+1, 2], dtype=np.float64)

    for it in range(int(time/step)+1):
        IN_TDDFT_TIME[it, 0] = it*step
        IN_TDDFT_TIME[it, 1] = it * step * strength

    if save_:
        np.savetxt('IN.TDDFT_TIME', IN_TDDFT_TIME)

    return IN_TDDFT_TIME

#--------------------------------------------------------------------
'''
    t<time1: ftddft=a*time^2
    t>time1: ftddft=strength*time
'''
def Gen_QuadraticField(time1, time, step, strength, save_=True):

    a = strength/(2*time1)

    IN_TDDFT_TIME = np.zeros([int(time/step)+1, 2], dtype=np.float64)
   
    for it1 in range(int(time1/step)+1):
        IN_TDDFT_TIME[it1, 0] = it1*step
        ftddft                = a * (it1*step)**2
        IN_TDDFT_TIME[it1, 1] = ftddft

    b = a*time1**2 - strength * time1 

    for it2 in range(int(time1/step)+1, int(time/step)+1):
        IN_TDDFT_TIME[it2, 0] = it2*step
        ftddft                = strength * it2*step + b
        IN_TDDFT_TIME[it2, 1] = ftddft 

    if save_:
        np.savetxt('IN.TDDFT_TIME', IN_TDDFT_TIME)

    return IN_TDDFT_TIME

#--------------------------------------------------------------------
def Plot(IN_TDDFT_TIME):

    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)

    time   = IN_TDDFT_TIME[:,0]
    ftddft = IN_TDDFT_TIME[:,1]

    plt.plot(time, ftddft, lw=2)

    plt.axvline(x=500, ymin=0, ymax=np.max(ftddft), color='red')

    plt.show()    


#--------------------------------------------------------------------
def main():

    gen_mode = 'linear' # linear; quadratic_linear

    if gen_mode == "linear":
        time     = 1000 # fs
        step     = 0.1  # fs
        # 1 MV/cm = 0.01 V/Ang = 10 mV/Ang
        # 1 kV/cm = 1*10^-5 V/Ang
        strength = 0.1 # strength*|A|(0.1) eV/Ang or V/Ang
        IN_TDDFT_TIME = Gen_LinearField(time, step, strength)
        Plot(IN_TDDFT_TIME)

    if gen_mode == 'quadratic_linear':
        time1     = 500  # fs
        time      = 1000 # fs
        step      = 0.1  # fs
        # 1 MV/cm = 0.01 V/Ang = 10 mV/Ang
        # 1 kV/cm = 1*10^-5 V/Ang
        strength  = 0.1 # strength*|A|(0.1) eV/Ang or V/Ang        
        IN_TDDFT_TIME = Gen_QuadraticField(time1, time, step, strength)
        Plot(IN_TDDFT_TIME)

if __name__ == "__main__":
    main()
