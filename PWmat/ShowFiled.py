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
    show the electronic filed for pwamt
    TDDFT_TIME
    itype_time = 1 or 2
    check the manual of Eneglish version
'''

import numpy as np
import scipy.integrate as si
import matplotlib as mpl
import matplotlib.pyplot as plt

#-----------------------------------------------------------
def functions(x):

    global b_param

    b1, b2, b3, b4, b5 = b_param[:]
    
    return b1 * np.exp(-((x-b2)**2 / (b3)**2)) * np.sin(b4*x/np.pi + b5)

#-----------------------------------------------------------
'''
    f = b1\cdot e^{-\frac{\left(t-b2\right)^{2}}{b3^{2}}}\cdot\sin\left(b4\cdot t+b5\right)
'''
def itype_time1(time):

    global b_param

    b1, b2, b3, b4, b5 = b_param[:]

    ntime = np.shape(time)[0]
    f     = np.zeros([ntime])

    for i in range(ntime):
        it   = time[i]
        f[i] = b1 * np.exp(-((it-b2)**2 / (b3)**2)) * np.sin(b4*it/np.pi + b5)
    
    return f

#----------------------------------------------------------
'''
    f = \int_{0}^{t}b1\cdot e^{-\frac{\left(t-b2\right)^{2}}{b3^{2}}}\cdot\sin\left(b4\cdot t+b5\right)dt
'''
def itype_time2(time):

    ntime = np.shape(time)[0]
    f     = np.zeros([ntime])

    for i in range(ntime):
        it   = time[i]
        f[i] = si.quad(functions, 0, it)[0]

    return f

def Plot(time, f):

    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)

    plt.plot(time, f, lw=2)

    plt.show()


#-----------------------------------------------------------
def main():

    itype_time = 2 # 1 or 2 in TDDFT_TIME
    # [b1, b2, b3, b4, b4] 
    global b_param
    #b_param    = [0.18, 50, 25, 5.027, 1.5708]
    b_param    = [1, 1, 100000, 0, np.pi/2]

    # get time, 0.1fs
    time       = np.array([i*0.1 for i in range(10000)])

    # get f
    if itype_time == 1:
        f = itype_time1(time)
    elif itype_time == 2:
        f = itype_time2(time)
    else:
        exit('Error itype_time! 1 or 2')

    # plot
    Plot(time, f)


if __name__ == "__main__":
    main()

