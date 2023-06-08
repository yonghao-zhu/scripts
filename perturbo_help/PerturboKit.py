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

# a perturbo code tool

import numpy as np
import os
np.set_printoptions(precision=15)
# task 3
import matplotlib as mpl
import matplotlib.pyplot as plt

# task 4
from xml.dom.minidom import parse

# task 5
import linecache

############################################################################
# for phonon calculation
## task 1
calphonon = {}
calphonon['TorF']      = False # True or False
calphonon['ph-x']      = 47 # start 1
# prefix = ./tmp/
calphonon['ph-ref.in'] = 'ph-ref.in'

############################################################################
# for nscf calculation
## task 2
kptgenerater = {}
kptgenerater['TorF']    = False # True or False
kptgenerater['mp_grid'] = [8, 8, 8]
# kpt.txt and kpt_weight.txt
kptgenerater['_save']   = True # True or False

############################################################################
# for phonon dispersion plot
# from cm^-1 to meV (1240*cm^-1 = 10^7eV = 10^4meV)
## task 3
phdisp = {}
phdisp['TorF']      = False # True or False
phdisp['freq']      = 'graphene.freq.gp'
phdisp['_plot']     = False 
phdisp['ene_range'] = [0, 250, 50] # meV
phdisp['high_kpt']  = ['Γ', 'M', 'K', 'Γ']
phdisp['insert_k']  = 50
phdisp['_save']     = False

############################################################################
# transfrom dyn.xml to dyn
## task 4
xml2dat_p = {}
xml2dat_p['TorF']      = False # True or False
xml2dat_p['num_files'] = 21
xml2dat_p['prefix']    = "graphene"

############################################################################
# read and save ShengBTE log file
## task 5
GetSheng = {}
GetSheng['TorF']     = True # True or False
GetSheng['nfiles'] = 56 # ShengBTE outputs

############################################################################
#----------------------------- functions -----------------------------------
############################################################################
# task 1
def CalPhonon(calphonon):

	'''
	ph-ref.in:
	Phonons on a uniform grid
	&inputph
	  verbosity='debug'
	  tr2_ph=1.0d-17
	  prefix='graphene'
	  ldisp=.true.
	  lqdir = .true.
	  outdir='./tmp'
	  fildyn  = 'graphene.dyn.xml'
	  fildvscf = 'dvscf'
	  nq1=18, nq2=18, nq3=1,
	  !nk1=60, nk2=60, nk3=1
	  start_q = 1
	  last_q = 1
	/
	'''

	ph_ref_in = calphonon['ph-ref.in']

	with open(ph_ref_in, 'r') as f:
		ph_ref_in_ = f.readlines()

	for i in range(1, calphonon['ph-x']+1):
		# mkdir ph-x
		if not os.path.isdir('ph-%s' %i):
			print('mkdir ps-%s...' %i)
			os.system('mkdir ph-%s' %i)

		# rewrite ph-ref.in
		with open('./ph-%s/ph.in' %i, 'w+') as f:
			for lines in ph_ref_in_:
				if 'outdir' in lines:
					f.writelines('  outdir = \'./tmp\'' + '\n')
				elif 'start_q' in lines:
					f.writelines('  start_q = %s' %i + '\n')
				elif 'last_q' in lines:
					f.writelines('  last_q = %s' %i + '\n')
				else:
					f.writelines(lines)

if calphonon['TorF']:
	if not os.path.isfile(calphonon['ph-ref.in']):
		print('No %s!' %calphonon['ph-ref.in'])
	else:
		CalPhonon(calphonon)

############################################################################
# task 2
def KptGenerater(kptgenerater):
	mp_grid = kptgenerater['mp_grid']

	delta_x = 1/mp_grid[0]
	delta_y = 1/mp_grid[1]

	# 2D
	if mp_grid[2] == 1:
		delta_z = 0
	# 3D
	else:
		delta_z = 1/mp_grid[2]
		print(delta_z)

	KPT = []
	KPT_weight = []

	weight = 1/(mp_grid[0]*mp_grid[1]*mp_grid[2])

	for x in range(0, mp_grid[0]):
		for y in range(0, mp_grid[1]):
			# 2D
			if mp_grid[2] == 1:
				KPT.append([x * delta_x, y * delta_y, delta_z])
				KPT_weight.append([x * delta_x, y * delta_y, delta_z, weight])
			# 3D
			else:
				for z in range(0, mp_grid[2]):
					KPT.append([x * delta_x, y * delta_y, z * delta_z])
					KPT_weight.append([x * delta_x, y * delta_y, z * delta_z, weight])

	KPT = np.array(KPT)
	KPT_weight = np.array(KPT_weight)

	print('KPT.shape=', KPT.shape)

	if kptgenerater['_save']:
		np.savetxt('kpt.txt', KPT, fmt='%0.8f')
		np.savetxt('kpt_weight.txt', KPT_weight, fmt='%0.8f')

if kptgenerater['TorF']:
	KptGenerater(kptgenerater)

############################################################################
# task 3
def PhDisp(phdisp):
	'''
		k_point mode1 (cm^-1) mode2 ...
	'''

	# read phdisp['freq']
	freq = np.loadtxt(phdisp['freq'])
	print('%s.shape:' %phdisp['freq'], freq.shape, '(kpt, modes+1)')

	# plot phonon dispersion
	if phdisp['_plot'] or phdisp['_save']:

		plt.rc('font',family='Times New Roman')
		mpl.rcParams['xtick.direction'] = 'in'
		mpl.rcParams['ytick.direction'] = 'in'

		x_kpt = freq[:, 0]
		# cm^-1 --> meV
		y_freq = freq[:, 1:] * 1240 / 10000

		plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)

		# color = black
		#plt.plot(x_kpt, y_freq, c='black', lw=3)
		# none color
		plt.plot(x_kpt, y_freq, lw=4)

		npath = len(phdisp['high_kpt']) - 1
		high_kpt_index = [phdisp['insert_k']*i for i in range(npath)]
		
		y_range = phdisp['ene_range']

		plt.xlim([0,max(x_kpt)]); plt.ylim(y_range[:2])
		plt.ylabel("Phonon energy (meV)",size=28)

		x_h = [x_kpt[i] for i in high_kpt_index]

		for i in x_h:
			plt.axvline(x=i,c="black",lw=1.5)

		x_h.append(x_kpt[-1])
		plt.xticks(x_h, phdisp['high_kpt'], size=20)

		length = int((max(y_range[:2]) - min(y_range[:2])) / y_range[2])

		plt.yticks([min(y_range[:2]) + i*y_range[2] for i in range(length+1)],
		    [min(y_range[:2]) + i*y_range[2] for i in range(length+1)],size=20)

		if phdisp['_save']:
			plt.savefig('ph_dispersion.png', dpi=600, format='png')

			meV_data = np.zeros_like(freq)

			meV_data[:, 0] = freq[:, 0]
			meV_data[:, 1:] = freq[:, 1:] * 1240 / 10000

			np.savetxt('phonon_dispresion.txt', meV_data, fmt='%.6f')

		if phdisp['_plot']:
			plt.show()

if phdisp['TorF']:
	if os.path.isfile(phdisp['freq']):
		print('Running PhDisp......')
		PhDisp(phdisp)
	else:
		print('No %s! Exiting......' %phdisp['freq'])

############################################################################
# task 4
'''
    s = "1.000000000  0.000000570  0.000000000"
'''
def out(s):
    if len(s.split()) <= 1:
        exit("Error out!")

    new = ""
    for i in s.split():
        if "-" not in i:
            new += "  " + i
        else:
            new += " " + i

    return new

#-----------------------------------------------------------
'''
    read header from dyn1.xml
    return header
'''
def readheader(filename):
    '''
    Dynamical matrix file
    Phonons calculation                                                        
    1    2   0   4.6593243   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
    Basis vectors
      1.000000000    0.000000570    0.000000000
     -0.499999506    0.866025688    0.000000000
      0.000000000    0.000000000    6.083691499
           1  'C   '    10947.0833707051     
    1    1     -0.0000447471      0.0000258347      2.7078472390
    2    1      0.0000400691      0.5773277217      2.7078550432

     Dynamical  Matrix in cartesian axes
    '''

    massfactor = 1.8218779*6.022E-4

    header = ['Dynamical matrix file', 'Phonons calculation ']

    fcxml = parse(filename)
    data  = fcxml.documentElement

    # get natoms
    line   = data.getElementsByTagName('NUMBER_OF_ATOMS')[0].childNodes[0].nodeValue.split()
    natoms = int(line[0])

    # get ntype
    line   = data.getElementsByTagName('NUMBER_OF_TYPES')[0].childNodes[0].nodeValue.split()
    ntypes = int(line[0])    

    # get CELL_DIMENSIONS
    line   = data.getElementsByTagName('CELL_DIMENSIONS')[0].childNodes[0].nodeValue.split()
    cell   = [float(i) for i in line]
    cell   = "  %s  %s  0  %0.7f  %0.7f  %0.7f  %0.7f  %0.7f  %0.7f" %(ntypes, natoms, \
                                                                    cell[0], cell[1], cell[2], \
                                                                    cell[3], cell[4], cell[5])

    header.append(cell); header.append("Basis vectors")

    # get AT
    line   = data.getElementsByTagName('AT')[0].childNodes[0].nodeValue.split()
    cell   = [float(i) for i in line]
    cell_1 = "     " + out("%10.9f  %10.9f  %10.9f" %(cell[0], cell[1], cell[2]))
    cell_2 = "     " + out("%10.9f  %10.9f  %10.9f" %(cell[3], cell[4], cell[5]))
    cell_3 = "     " + out("%10.9f  %10.9f  %10.9f" %(cell[6], cell[7], cell[8]))

    header.append(cell_1); header.append(cell_2); header.append(cell_3)

    # masses
    for i in range(1, ntypes+1):
        type_name = 'TYPE_NAME.%s' %i
        mass_name = 'MASS.%s' %i
        line_type = data.getElementsByTagName(type_name)[0].childNodes[0].nodeValue.split()[0]
        line_mass = data.getElementsByTagName(mass_name)[0].childNodes[0].nodeValue.split()[0]
        mass      = str('%.10f' % (float(line_mass)/massfactor))
        header.append("           %s " %i + " '" + line_type + "'  " + mass)

    for i in range(1, natoms+1):
        name = "ATOM.%s" %i
        line = data.getElementsByTagName(name)[0].attributes.items()
        type = int(line[1][1])
        tau  = line[2][1].split()
        tau  = "%0.10f  %0.10f  %0.10f" %(float(tau[0]), float(tau[1]), float(tau[2]))
        tmp  = "    " + out("%s  " %i + "%s  " %type + tau)

        header.append(tmp)

    return header

#-----------------------------------------------------------
'''
    read dynamical mat
    return dynamat
'''
def readdyna(filename):

    fcxml = parse(filename)
    data  = fcxml.documentElement

    dynmat = []

    # get natoms
    line   = data.getElementsByTagName('NUMBER_OF_ATOMS')[0].childNodes[0].nodeValue.split()
    natoms = int(line[0])

    # get NUMBER_OF_Q
    line  = data.getElementsByTagName('NUMBER_OF_Q')[0].childNodes[0].nodeValue.split()
    num_q = int(line[0])

    # get dynamical mat
    for i in range(1, 1+num_q):
        tmp       = ["     Dynamical  Matrix in cartesian axes", "\n"]
        dyna_name = "DYNAMICAL_MAT_.%s" %i
        dyna_line = data.getElementsByTagName(dyna_name)[0].getElementsByTagName("Q_POINT")
        qpoint    = dyna_line[0].childNodes[0].nodeValue.split()
        qpoint    = "     q = (   %0.10f   %0.10f   %0.10f  )" %(float(qpoint[0]), \
                                                            float(qpoint[1]), \
                                                            float(qpoint[2]))
        tmp.append(qpoint); tmp.append("\n")

        for iatom in range(1, 1+natoms):
            for jatom in range(1, 1+natoms):
                pair      = "    %s    %s" %(iatom, jatom)
                pair_name = "PHI.%s.%s" %(iatom, jatom)
                dyna_mat  = data.getElementsByTagName(dyna_name)[0].getElementsByTagName(pair_name)
                dyna_mat  = dyna_mat[0].childNodes[0].nodeValue.split('\n')
                if ',' in dyna_mat[0]:
                    dyna_mat  = [i.split(',') for i in dyna_mat]
                else:
                    dyna_mat  = [i.split() for i in dyna_mat]
                del dyna_mat[0]; del dyna_mat[-1]

                tmp.append(pair)
                for m in range(3):
                    dyna_1 = out("%0.8f  %0.8f" %(float(dyna_mat[3*m+0][0]), float(dyna_mat[3*m+0][1])))
                    dyna_2 = out("%0.8f  %0.8f" %(float(dyna_mat[3*m+1][0]), float(dyna_mat[3*m+1][1])))
                    dyna_3 = out("%0.8f  %0.8f" %(float(dyna_mat[3*m+2][0]), float(dyna_mat[3*m+2][1])))
                    dyna   = dyna_1 + "  " + dyna_2 + "  " + dyna_3
                    tmp.append(dyna)
        dynmat.append(tmp)
    return dynmat

#-----------------------------------------------------------
'''
    read dynamical omega
    return omega
'''
def readomega(filename):

    OMEGA = []

    fcxml = parse(filename)
    data  = fcxml.documentElement

    # get natoms
    line   = data.getElementsByTagName('NUMBER_OF_ATOMS')[0].childNodes[0].nodeValue.split()
    natoms = int(line[0])

    # get omega
    for i in range(1, 1+3*natoms):
        tmp               = []
        omega_name        = "OMEGA.%s" %i
        displacement_name = "DISPLACEMENT.%i" %i
        
        omega = data.getElementsByTagName(omega_name)[0].childNodes[0].nodeValue.split()
        omega = "     freq (    %s) =      %0.6f [THz] =       %0.6f [cm-1]" %(i, float(omega[0]), float(omega[1]))
        tmp.append(omega)

        displ  = data.getElementsByTagName(displacement_name)[0].childNodes[0].nodeValue.split('\n')
        del displ[0]; del displ[-1]

        if ',' in displ[0]:
            a1     = float(displ[0].split(',')[0]); a2     = float(displ[0].split(',')[1])
            a3     = float(displ[1].split(',')[0]); a4     = float(displ[1].split(',')[1])
            a5     = float(displ[2].split(',')[0]); a6     = float(displ[2].split(',')[1])

            b1     = float(displ[3].split(',')[0]); b2     = float(displ[3].split(',')[1])
            b3     = float(displ[4].split(',')[0]); b4     = float(displ[4].split(',')[1])
            b5     = float(displ[5].split(',')[0]); b6     = float(displ[5].split(',')[1])
        else:
            a1     = float(displ[0].split()[0]); a2     = float(displ[0].split()[1])
            a3     = float(displ[1].split()[0]); a4     = float(displ[1].split()[1])
            a5     = float(displ[2].split()[0]); a6     = float(displ[2].split()[1])

            b1     = float(displ[3].split()[0]); b2     = float(displ[3].split()[1])
            b3     = float(displ[4].split()[0]); b4     = float(displ[4].split()[1])
            b5     = float(displ[5].split()[0]); b6     = float(displ[5].split()[1])            

        displ1 = " ( %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f )" %(a1, a2, a3, a4, a5, a6)
        displ2 = " ( %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f )" %(b1, b2, b3, b4, b5, b6)
        
        tmp.append(displ1); tmp.append(displ2)

        OMEGA.append(tmp)
    return OMEGA

#-----------------------------------------------------------
'''
    write dynamical data
'''
def writedyn(num_files, prefix):

    for i in range(1, num_files+1):
        filename_xml = "%s.dyn%i.xml" %(prefix, i)

        print(filename_xml)

        header  = readheader(filename_xml)
        dynamat = readdyna(filename_xml)
        omega   = readomega(filename_xml)

        # write files
        filename_dat = "./dyn/%s.dyn%i" %(prefix, i)
        with open(filename_dat, 'w+') as f:
            # write header
            for i in header:
                f.writelines(i + '\n')
            # write dynamical mat
            for i in dynamat:
                f.writelines("\n")
                for j in i:
                    if j != "\n":
                        f.writelines(j + "\n")
                    else:
                        f.writelines(j)
            f.writelines("\n" + "     Diagonalizing the dynamical matrix" + "\n \n")
            f.writelines("     " + dynamat[0][2] + "\n \n")
            f.writelines(" **************************************************************************\n")
            # write omega
            for i in omega:
                for j in i:
                    f.writelines(j + "\n")
            f.writelines(" **************************************************************************\n")

if xml2dat_p['TorF']:

	writedyn(xml2dat_p["num_files"], xml2dat_p["prefix"])

############################################################################
# task 5: read and save ShengBTE log
'''
	read lines
	return Vp2, ijk(3), q_m(3)
'''
def ReadEle(filename, num):

	line    = linecache.getline(filename, num).split()
	myid    = int(line[1])
	i, j, k = int(line[2]), int(line[3]), int(line[4])
	ijk     = np.array([i, j, k])

	q       = int(line[5])
	qprime  = int(line[6])
	qdprime = int(line[7])
	q_m     = np.array([q, qprime, qdprime])  

	Vp2          = float(line[8])

	return Vp2, ijk, q_m

#--------------------------------------------------------------------
'''
	N_plus : number of allowed absorption processes
	N_minus: number of allowed emission processes.
	maxsize = max(N_plus, N_minus)
	NList  : reduced qpoints
	Nbands : 3*natoms
'''
def readVp(nfiles):
	
	with open('log.dat', 'r') as f:
		for line in f.readlines():
			# get Nbands
			if 'Nbands' in line:
				Nbands = int(line.split()[-1]) 
			# get NList
			if 'NList' in line:
				NList = int(line.split()[-1])
			# get maxsize
			if 'maxsize' in line:
				maxsize = int(line.split()[-1])

	# check 
	print('NList    = ', NList)
	print('maxsize  = ', maxsize)
	print('Nbands   = ', Nbands)

	Vp_plus_file   = open('Vp3_plus.txt', 'w+')
	Vp_minus_file  = open('Vp3_minus.txt', 'w+')
	ijk_plus_file  = open('ijk3_plus.txt', 'w+')
	ijk_minus_file = open('ijk3_minus.txt', 'w+')
	q_plus_file    = open('q3_plus.txt', 'w+')
	q_minus_file   = open('q3_minus.txt', 'w+')

	Vp_plus_file.writelines(str(NList) + '  ' + str(Nbands) + '  ' + \
			         str(maxsize) + '\n')
	Vp_minus_file.writelines(str(NList) + '  ' + str(Nbands) + '  ' + \
			         str(maxsize) + '\n')

	plus_num = 0; minus_num = 0
	n_tot = 0
	for i in range(nfiles):
		filename = "./out/%d.txt" %(i+1)
		file = open(filename, 'r')
		for num, line in enumerate(file):
			# mm, N_plus_count, i, j, k, list(ll), ii, ss, 
			# real(dconjg(Vp_plus_matrix_reduce(mm, N_plus_count))*Vp_plus_matrix_reduce(mm, N_plus_count))
			# Vp_plus_matrix
			if 'Vp_m:' in line:
				Vp2, ijk, q_m = ReadEle(filename=filename, num=num+1)
				Vp_minus_file.writelines('%0.16E \n' %Vp2)
				ijk_minus_file.writelines(str(ijk[0]) + '  ' + str(ijk[1]) + '  ' + str(ijk[2]) + '\n')
				q_minus_file.writelines(str(q_m[0]) + '  ' + str(q_m[1]) + '  ' + str(q_m[2]) + '\n')
				minus_num += 1
				n_tot += 1

			# Vp_minus_matrix
			if 'Vp_p:' in line:
				Vp2, ijk, q_m = ReadEle(filename=filename, num=num+1)
				Vp_plus_file.writelines('%0.16E \n' %Vp2)
				ijk_plus_file.writelines(str(ijk[0]) + '  ' + str(ijk[1]) + '  ' + str(ijk[2]) + '\n')
				q_plus_file.writelines(str(q_m[0]) + '  ' + str(q_m[1]) + '  ' + str(q_m[2]) + '\n')
				plus_num += 1
				n_tot += 1

			if n_tot % 100000 == 0:
				print(n_tot)

		file.close()

	Vp_plus_file.close()
	Vp_minus_file.close()
	ijk_plus_file.close()
	ijk_minus_file.close()
	q_plus_file.close()
	q_minus_file.close()

	print("plus_num  = ", plus_num)
	print("minus_num = ", minus_num)

if GetSheng['TorF']:
    readVp(nfiles=GetSheng['nfiles'])

############################################################################



############################################################################