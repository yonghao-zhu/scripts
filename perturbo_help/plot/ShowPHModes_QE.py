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
    show phonon modes
    open the xsf with VESTA

    inputs: a POSCAR with vasp format and matdyn.modes with qe format

    email : yonghao_zhu@163.com
    github: https://github.com/Crazy-Rookie and https://github.com/yonghao-zhu
    
    2022-10-05
'''

import numpy as np
import os
import linecache

#-----------------------------------------------------------------------------
'''
    get positions from POSCAR with vasp format

    return latt, pos
'''
def GetPositions(poscar):

    # get latt
    latt = np.zeros([3, 3], dtype=np.float64)
    for i in range(3, 6):
        line = linecache.getline(poscar, i).split()
        for j in range(3):
            latt[i-3, j] = float(line[j])
    
    # get natoms
    num_atom = [int(i) for i in linecache.getline(poscar, 7).split()]
    num_type = linecache.getline(poscar, 6).split()
    natoms   = sum(num_atom)

    types = []
    for i in range(len(num_type)):
        types += [num_type[i]] * num_atom[i]

    print(" --> The number of atoms: ", natoms)

    print(num_atom)
    print(types)

    # get pos
    pos = np.zeros([natoms, 3], dtype=np.float64)
    for i in range(9, 9+natoms):
        line = linecache.getline(poscar, i).split()
        for j in range(3):
            pos[i-9, j] = float(line[j])
    
    if "D" in linecache.getline(poscar, 8):
        pos = np.dot(pos, latt)

    return latt, pos, types

#-----------------------------------------------------------------------------
'''
    get vibrations vectors from matdyn.modes with qe format
    iq = 0 for gamma point

    u + ui
    PHonon/FD/fd_ef.f90

    return vectors
'''
def GetVec(matdyn, natoms, iq_=0):

    nmodes = natoms * 3

    # get the number of qpoints
    qpts = 0
    with open(matdyn, 'r') as f:
        for i in f.readlines():
            if "q =" in i:
                qpts += 1
    
    print(" --> The number of qpoints: ", qpts)

    vectors = np.zeros([qpts, nmodes, natoms, 3], dtype=np.float64)

    for iq in range(qpts):
        lines = 3 + iq * ( (natoms+1) * nmodes + 5 )
        for imode in range(nmodes):
            for iatom in range(natoms):
                tmp  = 3 + lines + imode * (natoms+1) + iatom
                disp = linecache.getline(matdyn, tmp).split()
                disp = [float(disp[i]) for i in range(1, 7)]
                real = np.array([disp[0], disp[2], disp[4]])
                imag = np.array([disp[1], disp[3], disp[5]])
                
                vectors[iq, imode, iatom] = real

    return vectors[iq_]

#-----------------------------------------------------------------------------
'''
    save .xsf files
'''
def GetXsf(matdyn, poscar, m_select):

    latt, pos, types = GetPositions(poscar)
    
    natoms  = pos.shape[0]
    nmodes  = 3*natoms
    vectors = GetVec(matdyn, natoms, iq_=0)

    header  = ["CRYSTAL", "PRIMVEC", "PRIMCOORD", "%d  1" % natoms]
    
    lattice = []
    for i in range(3):
        tmp = "  %0.10f    %0.10f    %0.10f" %(latt[i, 0], latt[i, 1], latt[i, 2])
        lattice.append(tmp)

    for imode in range(len(m_select)):
        vec      = vectors[m_select[imode] - 1]
        filename = "gamma_%s.xsf" %m_select[imode]

        print(" --> filename: ", filename)

        with open(filename, "w+") as f: 

            f.writelines("CRYSTAL \n" + "PRIMVEC \n" + lattice[0] + "\n" + \
                                                       lattice[1] + "\n" + \
                                                       lattice[2] + "\n")
            f.writelines("PRIMCOORD \n" + "%d  1\n" % natoms)

            for iatom in range(natoms):
                tmp = types[iatom] + "  " + "%0.10f  %0.10f  %0.10f" %(pos[iatom][0],
                                                                       pos[iatom][1],
                                                                       pos[iatom][2])
                tmp += "  %0.6f  %0.6f  %0.6f" %(vec[iatom, 0],
                                                 vec[iatom, 1],
                                                 vec[iatom, 2])
                
                f.writelines(tmp + "\n")

#-----------------------------------------------------------------------------
def main():

    matdyn   = "../matdyn.modes"
    poscar   = "./POSCAR.vasp"
    # selected modes, starts 1
    m_select = [i for i in range(1, 10)]

    GetXsf(matdyn, poscar, m_select) 

if __name__ == "__main__":
    main()
