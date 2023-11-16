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
def CheckFile(filename):

    if type(filename) == str:
        if not os.path.isfile(filename):
            exit('No %s File! Exiting...' %filename)
        else:
            print('-->%s is ready!' %filename)
    elif type(filename) == list:
        for iname in filename:
            if not os.path.isfile(iname):
                exit('No %s File! Exiting...' %filename)
            else:
                print('-->%s is ready!' %iname)
    else:
        exit('I do not know the type of %s! Exiting...' %filename)               

#--------------------------------------------------------------------
def PlotOccEne(plot_file, chose_orbital, legend_, color_):

    # check file
    CheckFile(filename=plot_file)

    # read file
    occene = np.loadtxt(plot_file)
    print('--> occene.shape = ', occene.shape, '(times/fs, norbitals+1)')

    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.18)

    for i in range(len(chose_orbital)):
        if len(chose_orbital[i]) == 1:
            plt.plot(occene[:,0]/1000, occene[:, np.array(chose_orbital[i])], \
                     label=legend_[i], color=color_[i], lw=2)
        elif len(chose_orbital[i]) > 1:
            plt.plot(occene[:,0]/1000, occene[:, np.array(chose_orbital[i][0])], \
                     label=legend_[i], color=color_[i], lw=2)
            for j in range(1, len(chose_orbital[i])):
                plt.plot(occene[:,0]/1000, occene[:, np.array(chose_orbital[i][j])], \
                         color=color_[i], lw=2)
        else:
            exit('Error chose_orbital! Exiting...')           

    plt.legend(loc='upper right', prop={'size': 14})

    if plot_file == 'plot.TDDFT.DOS':
        plt.xlabel('Time (ps)', size=25); plt.ylabel('Occupation', size=25)
        plt.xlim(0, 1); plt.ylim(-0.2, 2.2)
    else:
        plt.xlabel('Time (ps)', size=25); plt.ylabel('Energy (eV)', size=25)
        plt.xlim(0, 1); plt.ylim(np.min(occene[:,1:])-0.5, np.max(occene[:,1:])+0.5)
    
    plt.xticks(fontsize=20); plt.yticks(fontsize=20)

    plt.show()    

#--------------------------------------------------------------------
def PlotRMSD(rmsd_file, legend_, color_):

    # check file
    CheckFile(rmsd_file)

    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.18)

    max_ = 0
    for i in range(len(rmsd_file)):
        rmsd_i = np.loadtxt(rmsd_file[i], comments=['frame', '0      NA'])
        if np.max(rmsd_i[:, 1]) > max_:
            max_ = np.max(rmsd_i[:, 1])

        plt.plot(rmsd_i[:,0]/1000, rmsd_i[:,1], label=legend_[i], color=color_[i], lw=2)

    plt.legend(loc='upper right', prop={'size': 14})

    plt.xlabel('Time (ps)', size=25); plt.ylabel('RMSD (Å)', size=25)
    plt.xlim(0, 1); plt.ylim(0, max_+0.1)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20)

    plt.show()   

#--------------------------------------------------------------------
def PlotMDSTPES(mdsteps_f, param):

    # check file
    CheckFile(mdsteps_f)

    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.18)

    # plot
    time = np.genfromtxt(mdsteps_f)[:, 1].astype(float)

    if param == 'etot':
        etot = np.genfromtxt(mdsteps_f)[:, 4].astype(float)
        plt.plot(time/1000, etot-etot[0], lw=2)
        plt.ylabel('Total Energy (eV)', size=25)
        plt.ylim(-2, 50)

    elif param == 'temp':
        temp = np.genfromtxt(mdsteps_f)[:, 7].astype(float)
        plt.plot(time/1000, temp, lw=2)
        plt.ylabel('Temperature (K)', size=25)
        plt.ylim(0, 600)

    elif param == 'scf':
        scfs = np.genfromtxt(mdsteps_f)[:, -5].astype(int)
        plt.plot(time/1000, scfs, lw=2)
        plt.ylabel('SCF Steps', size=25)
        plt.ylim(1, max(scfs)+1)

    elif param == 'dE':
        dE = np.genfromtxt(mdsteps_f)[:, -9].astype(float)
        plt.plot(time/1000, dE, lw=2)
        plt.ylabel('dE', size=25)
        plt.ylim(min(dE), max(dE))

    elif param == 'dR':
        dRhoe = np.genfromtxt(mdsteps_f)[:, -7].astype(float)
        plt.plot(time/1000, dRhoe, lw=2)
        plt.ylabel('dRhoe', size=25)
        plt.ylim(min(dRhoe), max(dRhoe))

    else:
        exit('Exiting...')

    plt.xlabel('Time (ps)', size=25)
    plt.xlim(0, 1)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20)

    plt.show()    

#--------------------------------------------------------------------
def PlotTDDOS(itime, DosType, ene_rang, dos_rang):

    total_dos_f = './RunDos/%s/total_dos/DOS.totalspin_projected_ShiftFermi' %itime
    occ_dos_f   = './RunDos/%s/occ_dos/DOS.totalspin_projected_ShiftFermi' %itime

    CheckFile([total_dos_f, occ_dos_f])

    ene       = np.loadtxt(total_dos_f, comments=['#'])[:, 0]
    total_dos = np.loadtxt(total_dos_f, comments=['#'])[:, 1]
    occ_dos   = np.loadtxt(occ_dos_f, comments=['#'])[:, 1]

    if total_dos.shape != occ_dos.shape:
        exit('total_dos.shape != occ_dos.shape')

    unocc_dos = total_dos[:] - occ_dos[:]

    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.18)

    if DosType == 'total':
        plt.plot(ene, total_dos, lw=2)
    elif DosType == 'occ':
        plt.plot(ene, occ_dos, lw=2)
    elif DosType == 'unocc':
        plt.plot(ene, unocc_dos, lw=2)
    elif DosType == 'all':
        plt.plot(ene, unocc_dos, label='unocc', color='blue', lw=2)
        plt.plot(ene, occ_dos, label='occ', color='red', lw=2)
        #plt.plot(ene, total_dos, label='total', color='black', lw=2)

        plt.legend(loc='upper right', prop={'size': 14})

    plt.xlabel('Energy (eV)', size=25); plt.ylabel('DOS (eV^-1)', size=25)
    plt.xlim(ene_rang); plt.ylim(dos_rang)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20)

    plt.show() 

#--------------------------------------------------------------------
def ReadRHO(filename):

    CheckFile(filename)

    natm = int(linecache.getline(filename, 7).split()[0])
    grid = [int(i) for i in linecache.getline(filename, natm+11).split()]

    points = grid[0]*grid[1]*grid[2]
    if points % 6:
        nlines = int(points / 6)
    else:
        nlines = int(points / 6) + 1

    atom = []
    for i in range(1, natm + 17):
        atom.append(linecache.getline(filename, i))

    rho = []
    for i in range(nlines):
        rho_line = linecache.getline(filename, natm+17+i)
        for ii in rho_line.split():
            rho.append(float(ii))

    return atom, np.array(rho)

def WriteRHO(out_name, atom, rho_sum):

    with open(out_name, 'w+') as f:
        for ia in atom:
            f.writelines(ia)
        for i in range(0, rho_sum.shape[0], 6):
            irho = '    '.join(map(str, rho_sum[i:i+6]))
            f.writelines(irho+'\n')

def SumRHO(itimes, index_, save_=True):

    for it in itimes:
        path = './RunDos/%s/charge_density/' %it

        for i in range(len(index_)):
            out_name = path+'RHO.%s.xsf' %index_[i][0]
    
            if len(index_[i]) == 2:
                in_name = path+'RHO.%s.xsf' %index_[i][1]
                with open(out_name, 'w+') as f, open(in_name, 'r') as g:
                    f.writelines(g.readlines())
            else:
                atom, rho_sum = ReadRHO(filename=path+'RHO.%s.xsf' %index_[i][1])
                for j in index_[i][2:]:
                    atom, rho = ReadRHO(filename=path+'RHO.%s.xsf' %j)
                    rho_sum += rho
                WriteRHO(out_name, atom, rho_sum)

#--------------------------------------------------------------------


#--------------------------------------------------------------------
def main():

    params         = {}
    # SumRho;   PlotTDDos;   PlotOccEne
    # PlotRMSD; PlotMDSTPES; 
    params['mode'] = 'PlotOccEne'

    #----------------------------------------------------------------
    if params['mode'] == 'PlotOccEne':
        # plot.TDDFT.E
        # plot.TDDFT.DOS
        plot_file = './data/plot.TDDFT.DOS' 
        # start 1
        chose_orbital  = [
                          [1,2,3,4], [5,6,7], [8],[9,10,11,12,13,14,15,16]
                         ]
        legend_        = ['VB','VBM', 'CBM','CB']
        color_         = ['red','black','blue','orange']

        PlotOccEne(plot_file, chose_orbital, legend_, color_)

    #----------------------------------------------------------------
    if params['mode'] == 'PlotRMSD':
        # need RMSD
        rmsd_file = [
                     'RMSD-Ga.dat', 'RMSD-N.dat', 
                     'RMSD-All.dat'
                    ]
        legend_    = ['Ga', 'N', 'All']
        color_     = ['blue', 'red', 'black']
        PlotRMSD(rmsd_file, legend_, color_)

    #----------------------------------------------------------------
    if params['mode'] == 'PlotMDSTPES':
        # need MDSTEPS
        mdsteps_f = './data/MDSTEPS'
        # etot; temp; scf; dE; dR
        param     = 'etot'
        PlotMDSTPES(mdsteps_f, param)

    #----------------------------------------------------------------
    if params['mode'] == 'SumRho':
        # need xsf files
        itimes = ['dos-0010', 'dos-0020']
        index_ = [ 
                   ['RHO.vb.xsf', '1-1', '1-2', '1-3'],
                   ['RHO.def.xsf', '1-4', '1-5', '1-6'],
                   ['RHO.cb.xsf', '1-7']
                 ]

        SumRHO(itimes, index_, save_=True)

    #----------------------------------------------------------------
    if params['mode'] == 'PlotTDDos':
        itime    = 'dos-1000'
        DosType  = 'all' # total, occ, unocc, all
        ene_rang = [-3, 3]
        dos_rang = [0, 80]

        PlotTDDOS(itime, DosType, ene_rang, dos_rang)

if __name__ == '__main__':
    main()


