#! /usr/bin/python3
# -*- conding=UTF-8 -*-

'''
    1. I can correctly calculate distance between two atoms.
    2. Support only Direct type of POSCAR. 
    3. author: yonghao_zhu@163.com
'''

import numpy as np
import linecache, os

np.set_printoptions(precision=15)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def ReadPOSCAR(params):

    '''
    input files: POSCAR

    return pos (Direct), latt (shape=[(atoms, 3), (3, 3)])
    '''

    poscar = params["posfile"]

    atoms_kind = linecache.getline(poscar, 6).split()
    atoms_num  = [int(i) for i in linecache.getline(poscar, 7).split()]
    #test# print(atoms_kind,atoms_num)

    if params["print_"]:
        print('++++++++++++++++++++POSCAR+++++++++++++++++++++')
        print('',atoms_kind,'\n',atoms_num)

        for i in range(len(atoms_kind)):
            if i == 0:
                print('', atoms_kind[i], '-->', '1 -', atoms_num[i])

            if i > 0:
                print('', atoms_kind[i], '-->', sum(atoms_num[:i])+1, '-', 
                      sum(atoms_num[:i]) + atoms_num[i])
        print('++++++++++++++++++++++++++++++++++++++++++++++')

    comments = [linecache.getline(poscar, i).rstrip('\n') for i in range(1, 9)]
    pos_np   = np.loadtxt(poscar, comments=comments)
    latt     = np.array([[float(m) for m in linecache.getline(poscar, i).split()] 
                          for i in range(3, 6)])

    return pos_np, latt


def ModificationCoordinate(pos):

    '''
        Support Only Direct Type
        pos.shape = (3)
        
        return pos_new + supercell_ (shape=(27, 3))
    '''

    # move to the center box
    pos_new = np.zeros_like(pos)
    for m in range(3): 
        if pos[m] > 1:
            pos_new[m] = pos[m] - 1
        if pos[m] < 0:
            pos_new[m] = pos[m] + 1
        if pos[m] >= 0 and pos[m] <= 1:
            pos_new[m] = pos[m]

    # make supercell (3,3,3)
    supercell_ = []
    for m in range(-1,2):
        for n in range(-1,2):
            for q in range(-1,2):
                a = []
                a.append(m)
                a.append(n)
                a.append(q)
                supercell_.append(a)
    supercell_ = np.array(supercell_).reshape(27,3)

    return pos_new + supercell_


def Dimension(pos_Atom1, pos_Atom2, dimension):

    """
        pos_Atom1.shape = (3)
        pos_Atom2.shape = (27, 3)

        return pos_Atom1, pos_Atom2 (shape=[(3), (27, 3)])
    """

    # one dimension
    if dimension == "x":
        pos_Atom1[1:]    = 0
        pos_Atom2[:, 1:] = 0

    elif dimension == "y":
        pos_Atom1[0]    = 0
        pos_Atom1[2]    = 0
        pos_Atom2[:, 0] = 0
        pos_Atom2[:, 2] = 0

    elif dimension == "z":
        pos_Atom1[:2]    = 0
        pos_Atom2[:, :2] = 0

    # two dimension
    elif dimension == "xy":
        pos_Atom1[2]    = 0
        pos_Atom2[:, 2] = 0

    elif dimension == "xz":
        pos_Atom1[1]    = 0
        pos_Atom2[:, 1] = 0

    elif dimension == "yz":
        pos_Atom1[0]    = 0
        pos_Atom2[:, 0] = 0 

    return pos_Atom1, pos_Atom2


def Distance(pos, latt, dimension = "xyz", print_ = False):

    """
        input: pos.shape =(2, 3)
               latt.shape=(3, 3)

        return  distance (float)
    """
    nbonds   = pos.shape[0]
    distance = np.zeros(nbonds, dtype=float)

    for ib in range(nbonds):
        pos_Atom1 = pos[ib, 0].dot(latt)
        pos_Atom2 = ModificationCoordinate(pos=pos[ib, 1]).dot(latt)
        # dimension 
        pos_Atom1, pos_Atom2 = Dimension(pos_Atom1 = pos_Atom1, 
                                         pos_Atom2 = pos_Atom2,
                                         dimension = dimension)
        distance_     = np.linalg.norm(pos_Atom1 - pos_Atom2, axis=1)
        distance[ib]  = np.min(distance_)

        if print_:
            print('-->Atoms: start 1')
            print('-->Distance Atom1-Atom2 : %.04f' %distance[ib],'Ang')
            print('++++++++++++++++++++++++++++++++++++++++++++++')

    return distance


def ReadXDATCAR(xdatcar, TimeRange):

    """
        input: XDATCAR

        return distance (shape=(times))
    """

    # time range
    time_ini = TimeRange[0]
    time_fin = TimeRange[1]

    # read XDATCAR in time range
    lattice  = np.mat([[float(linecache.getline(xdatcar, i+1).split()[0]),
                        float(linecache.getline(xdatcar, i+1).split()[1]),
                        float(linecache.getline(xdatcar, i+1).split()[2])] for i in range(2,5)])

    comments = ['Direct configuration=']
    atoms    = int(sum([float(i) for i in linecache.getline(xdatcar, 7).split()]))

    for i in range(7):
        comments.append(linecache.getline(xdatcar, i+1))

    xdatcar = np.loadtxt(xdatcar, comments=comments).reshape(-1, atoms, 3)[time_ini-1: time_fin, :, :]

    return xdatcar, lattice, atoms


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def main():

    '''
    1. input files    : POSCAR(Direct), XDATCAR(Direct) with vasp format
    2. Atoms and Tiems: start 1
    '''

    params             = {}
    params['posfile']  = './XDATCAR' # POSCAR (0K) or XDATCAR (300K)
    params["dimension"]= 'xyz'     # x (1), y(2), z(3), xy, xz, yz, xyz
    params['AtomList'] = [
                           [75, 172],[75, 170],[75, 171],[75, 162],
                           [77, 172],[77, 173],[77, 174],[77, 164],
                           [85, 172],[85, 180],[85, 181],[85, 182],

                           [37, 141],[37, 131],[37, 133],[37, 164],
                           [45, 141],[45, 132],[45, 140],[45, 142],
                           [46, 141],[46, 143],[46, 151],[46, 174],

                           [52, 149],[52, 147],[52, 157],[52, 180],
                           [53, 149],[53, 140],[53, 148],[53, 150],
                           [54, 149],[54, 151],[54, 182],[54, 159],

                           [35, 139],[35, 129],[35, 131],[35, 162],
                           [43, 139],[43, 137],[43, 147],[43, 170],
                           [44, 139],[44, 130],[44, 138],[44, 140],
                         ]         # including, start 1
    params['TimeRange']= [1, 10000]    # only for md, including, start 1
    params["save_"]    = True      # only for md
    params["print_"]   = True      # True or False

    if 'XDATCAR' not in params['posfile']:
        pos_np, latt = ReadPOSCAR(params)

        # distance between two atoms
        nbonds = len(params['AtomList'])
        pos    = np.zeros([nbonds, 2, 3])
        for ib in range(nbonds):
            pos[ib, 0] = pos_np[params["AtomList"][ib][0]-1]
            pos[ib, 1] = pos_np[params["AtomList"][ib][1]-1]

        Distance(pos = pos, latt = latt, dimension = params["dimension"], print_=True)

    if 'XDATCAR' in params['posfile']:
        xdatcar, lattice, atoms = ReadXDATCAR(xdatcar   = params['posfile'], 
                                              TimeRange = params["TimeRange"])

        time_ini = params["TimeRange"][0]
        time_fin = params["TimeRange"][1]
        nbonds   = len(params['AtomList'])
        dis      = np.zeros([time_fin - time_ini + 1, 1+nbonds], dtype=float)

        # distance between A1 and A2
        for itime in range(time_fin - time_ini + 1):
            pos = np.zeros([nbonds, 2, 3])
            for ib in range(nbonds):
                pos[ib, 0] = xdatcar[itime, params['AtomList'][ib][0]-1, :]
                pos[ib, 1] = xdatcar[itime, params['AtomList'][ib][1]-1, :]

            dis[itime, 0]  = itime + 1
            dis[itime, 1:] = Distance(pos = pos, latt = lattice, dimension = params["dimension"])[:]

        if params["save_"]:
            np.savetxt("distance.txt", dis, fmt="%0.6f")

if __name__ == "__main__":
	main()