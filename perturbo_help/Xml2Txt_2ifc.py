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

import numpy as np
import os
from xml.dom.minidom import parse

'''
    update: 2023-04-19
    use this script to transform 2ifc.xml to 2ifc (txt format)
'''

xml_fname  = 'mos212121.fc.xml'
txt_fname  = 'mos212121.fc'
Leps_zstar = True # True or False

#-----------------------------------------------------------------------------
'''
    read xxx.ifc.xml with xml format

    return ph
'''
def ReadIFC2XML(ifc2_fname, Leps_zstar):

    assert os.path.isfile(ifc2_fname), exit('No %s! Exitting...' %ifc2_fname)
    print("===============* ReadIFC2XML *===============")

    fcxml = parse(ifc2_fname)
    
    # read nat
    tmp_line = fcxml.getElementsByTagName('NUMBER_OF_ATOMS')[0].childNodes[0].nodeValue.split()
    nat      = int(tmp_line[0])
    print('--1# info: natoms ', nat)

    # read ntype
    tmp_line = fcxml.getElementsByTagName('NUMBER_OF_TYPES')[0].childNodes[0].nodeValue.split()
    ntype    = int(tmp_line[0])
    print('--2# info: ntype ', ntype)

    # read celldm
    tmp_line = fcxml.getElementsByTagName('CELL_DIMENSIONS')[0].childNodes[0].nodeValue.split()
    celldm   = np.array([float(i) for i in tmp_line])

    # read lattice vector
    tmp_line = fcxml.getElementsByTagName('AT')[0].childNodes[0].nodeValue.split()
    at       = np.array([float(i) for i in tmp_line]).reshape([3, 3])
    print('--3# info: crystal lattice \n', at)

    # read mass
    amass = {}
    for i in range(1, ntype+1):
        type_name = 'TYPE_NAME.%s' %i
        mass_name = 'MASS.%s' %i
        line_type = fcxml.getElementsByTagName(type_name)[0].childNodes[0].nodeValue.split()[0]
        line_mass = fcxml.getElementsByTagName(mass_name)[0].childNodes[0].nodeValue.split()[0]
        amass[i]  = float(line_mass)
    print('--4# info: mass ', amass)

    # read tau
    tau  = np.zeros([nat, 3], dtype=np.float64)
    ityp = np.zeros([nat], dtype=np.int64)
    for i in range(1, nat+1):
        name      = "ATOM.%s" %i
        line      = fcxml.getElementsByTagName(name)[0].attributes.items()
        ityp[i-1] = int(line[1][1])
        tau[i-1]  = line[2][1].split()
    print('--5# info: tau, ityp \n', tau, '\n', ityp)

    # read omega, Bohr^3
    tmp_line = fcxml.getElementsByTagName('UNIT_CELL_VOLUME_AU')[0].childNodes[0].nodeValue.split()
    omega    = float(tmp_line[0])
    print('--3# info: omega = %0.4f' %omega)

    # read scell
    tmp_line = fcxml.getElementsByTagName('MESH_NQ1_NQ2_NQ3')[0].childNodes[0].nodeValue.split()
    scell    = np.array([int(tmp_line[0]), int(tmp_line[1]) ,int(tmp_line[2])])
    print('--6# info: scell ', scell)

    # read bg, reciprocal lattice
    tmp_line = fcxml.getElementsByTagName('BG')[0].childNodes[0].nodeValue.split()
    bg       = np.array([float(i) for i in tmp_line]).reshape([3, 3])
    print('--7# info: reciprocal lattice \n', bg)    

    if Leps_zstar:
        # read epsilon
        tmp_line = fcxml.getElementsByTagName('EPSILON')[0].childNodes[0].nodeValue.split()
        epsilon  = np.array([float(i) for i in tmp_line], dtype=np.float64).reshape([3, 3])
        print('--8# info: epsilon \n', epsilon)  
        # read zstar
        tmp_line = fcxml.getElementsByTagName('ZSTAR')[0].childNodes
        zstar    = np.zeros([nat, 3, 3], dtype=np.float64)
        for ia in range(nat):
            tmp_z     = tmp_line[2*ia+1].childNodes[0].nodeValue.split()
            zstar[ia] = np.array([float(i) for i in tmp_z]).reshape([3, 3])

    # read fc
    frc   = np.zeros([scell[0], scell[1], scell[2], 3, 3, nat, nat], dtype=np.float64)
    for na in range(nat):
        for nb in range(nat):

            for m3 in range(scell[2]):
                for m2 in range(scell[1]):
                    for m1 in range(scell[0]):

                        tag_name = 's_s1_m1_m2_m3.%d.%d.%d.%d.%d' %(na+1, nb+1, m1+1, m2+1, m3+1)
                        tmp_line = fcxml.getElementsByTagName(tag_name)[0]\
                                  .childNodes[1].childNodes[0].nodeValue.split()
                        ifc      = np.array([float(i) for i in tmp_line], dtype=np.float64).reshape([3, 3])
                        frc[m1, m2, m3, :, :, na, nb] = ifc.T

    # dict
    ph           = {}
    ph["ntype"]  = ntype
    ph["ifc"]    = frc
    ph["nat"]    = nat
    ph["celldm"] = celldm
    ph["at"]     = at
    ph["omega"]  = omega
    ph["amass"]  = amass
    ph["tau"]    = tau
    ph["ityp"]   = ityp
    ph["scell"]  = scell
    if Leps_zstar:
        ph['epsilon'] = epsilon
        ph['zstar']   = zstar

    return ph

#-----------------------------------------------------------------------------
'''
    NOTE: do not write epsilon and born effective charge!
'''
def WriteHeader(ph, txt_fname, Leps_zstar):

    AMU_SI           = 1.66053906660E-27 # Kg
    ELECTRONMASS_SI  = 9.1093837015E-31  # Kg
    AMU_AU           = AMU_SI / ELECTRONMASS_SI
    AMU_RY           = AMU_AU / 2.0

    ntype  = ph["ntype"]
    nat    = ph["nat"]
    celldm = ph["celldm"]
    at     = ph["at"] 
    amass  = ph["amass"]
    tau    = ph["tau"]  
    ityp   = ph["ityp"]   
    scell  = ph["scell"] 

    with open(txt_fname, 'w+') as f:
        line_1 = "  %d    %d   0  " %(ntype, nat)
        for i in celldm:
            line_1 += '%0.9f  ' %i
        f.writelines(line_1 + '\n')
        # lattice vector
        for i in range(3):
            line_2 = "%0.9f     %0.9f     %0.9f" %(at[i, 0],\
                                                   at[i, 1],\
                                                   at[i, 2])
            f.writelines(line_2.rjust(47) + '\n')
        # elements
        for i in range(ntype):
            line_3 = "%d  \'  \'    %0.10f" %(i+1, amass[i+1]*AMU_RY)
            f.writelines(line_3.rjust(38) + '\n')
        # tau
        for i in range(nat):
            line_4 = ("    %d    %d" %(i+1, ityp[i])).ljust(15) + \
                     ("%0.10f      %0.10f      %0.10f" %(tau[i, 0], tau[i, 1], tau[i, 2])).rjust(49)
            f.writelines(line_4.ljust(64) + '\n')
        # NOTE: F --> not change
        if Leps_zstar:
            f.writelines(' T \n')
            for i in range(3):
                eps    = ph['epsilon'] 
                line_6 = ('%0.7f' %eps[i, 0]).rjust(15) + ('%0.7f' %eps[i, 1]).rjust(15) +('%0.7f' %eps[i, 2]).rjust(15)  
                f.writelines(line_6+'\n')
            for i in range(nat):
                f.writelines('    %d' %(i+1)+'\n')
                for j in range(3):
                    zstar  = ph['zstar']
                    line_7 = ('%0.7f' %zstar[i, j, 0]).rjust(15) + ('%0.7f' %zstar[i, j, 1]).rjust(15) +\
                             ('%0.7f' %zstar[i, j, 2]).rjust(15)
                    f.writelines(line_7+'\n')
        else:
            f.writelines(' F\n')
        # scell
        line_5 = '%d  %d  %d' %(scell[0], scell[1], scell[2])
        f.writelines(line_5.rjust(12) + '\n')

#-----------------------------------------------------------------------------
'''
    write 2ifc
'''
def WriteIFC(ifc, nat, scell, txt_fname):

    assert os.path.isfile(txt_fname), exit('No %s! Exiting...' %txt_fname)

    #frc = np.zeros([scell[0], scell[1], scell[2], 3, 3, nat, nat], dtype=np.float64)
    with open(txt_fname, 'a+') as f:
        for ipol in range(3):
            for jpol in range(3):

                for na in range(nat):
                    for nb in range(nat):
                        tag = '%d   %d   %d   %d' %(ipol+1, jpol+1, na+1, nb+1)
                        f.writelines(tag.rjust(17)+'\n')

                        for n3 in range(scell[2]):
                            for n2 in range(scell[1]):
                                for n1 in range(scell[0]):        
                                    ifc_ = ifc[n1, n2, n3, ipol, jpol, na, nb]
                                    line = ('%d   %d   %d' %(n1+1, n2+1, n3+1)).rjust(12) +\
                                           ('%0.10E' %(ifc_)).rjust(20)
                                    f.writelines(line+'\n')

#-----------------------------------------------------------------------------
def Xml2Txt(xml_fname, txt_fname, Leps_zstar):

    assert os.path.isfile(xml_fname), exit('No %s! Exiting...' %xml_fname)

    ph = ReadIFC2XML(ifc2_fname = xml_fname, Leps_zstar=Leps_zstar)

    print("===============* Xml2Txt *===============")

    WriteHeader(ph=ph, txt_fname=txt_fname, Leps_zstar=Leps_zstar)

    WriteIFC(ifc=ph['ifc'], nat=ph['nat'], scell=ph['scell'], txt_fname=txt_fname)

Xml2Txt(xml_fname=xml_fname, txt_fname=txt_fname, Leps_zstar=Leps_zstar)
