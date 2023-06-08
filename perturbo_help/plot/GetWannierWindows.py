#! /usr/bin/python3
# -*- conding=UTF-8 -*-
#  .--,       .--,
# ( (  \.---./  ) )
#  '.__/o   o\__.'
#     {=  ^  =}
#      >  -  <
#     /  Zhu  \
#    //  Yong \\
#   //|  Hao  |\\
#   "'\       /'"_.-~^`'-.
#      \  _  /--'         `
#    ___)( )(___
#   (((__) (__)))  

'''
    to determine wannier energy windows

    yonghao_zhu@163.com
    2022/10/09 
    organization: https://github.com/Crazy-Rookie
    useful url:
    https://blog.sciencenet.cn/blog-2909108-1154273.html
    https://blog.sciencenet.cn/blog-2909108-1263724.html

    only for QE output files
    but easily extend to vasp format
'''

import os
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

def PlotBand(filename, band_index, erange, nband):

    assert os.path.isfile(filename), "No %s! Exiting..." %filename
    
    band = np.loadtxt(filename).reshape([nband, -1, 2])

    band  = band[np.array(band_index)]
    print("band.shape = ", band.shape)
    ene   = band[:, :, 1]
    kpath = band[:, :, 0]
    print(ene.shape)

    # print emin and emax
    for i in range(len(band_index)):
        emin = np.min(ene[i])
        emax = np.max(ene[i])
        print("-->band = ", band_index[i]+1, "  emin = ", emin, "  emax = ", emax)

    # plot
    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    plt.figure(figsize=(6,5))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)

    for i in range(len(band_index)):
        plt.plot(kpath[i], ene[i], c='black', zorder=0, lw=2)

    plt.xlim([min(kpath[0]), max(kpath[0])])
    plt.ylim([erange[0], erange[1]])

    plt.ylabel("Energy (eV)", size=20)
    plt.xticks([], size=15)

    length = int((max(erange[:2]) - min(erange[:2])) / erange[2])
    plt.yticks([min(erange[:2]) + i*erange[2] for i in range(length+1)],
        [min(erange[:2]) + i*erange[2] for i in range(length+1)],size=20)    

    plt.show()

#--------------------------------------------------------------------------
def main():

    filename   = "band.dat.gnu"
    erange     = [-5, 7, 2] # eV
    nband      = 120
    band_index = [i for i in range(64, 108)] # start 0

    PlotBand(filename=filename, band_index=band_index, erange=erange, nband=nband)

if __name__ == "__main__":
    main()

