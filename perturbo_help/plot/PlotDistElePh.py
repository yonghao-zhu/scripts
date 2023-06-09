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
	plot electron and phonon distributions along time evolution
	ShengBTE and Perturbo results are inputs.

	update: 2022-11-19
	author: yonghao_zhu@163.com
'''

# modules
import numpy as np
import h5py 
import matplotlib as mpl
import matplotlib.pyplot as plt 
import os

# inputs
params                = {}
params['mode']        = 'ele_dist' # ph_sheng, ph_dist, ele_dist, ele_band

# read from xxx_epwan.h5, lattice vector = alat * at
params['at']          = np.array([[1.0, 0.0, 0.0],
                                  [-0.5, 0.8660254, 0.0],
					              [0.0, 0.0, 4.87804878]])
# Bohr2Ang  = 5.29177208590000E-01					
params['alat']        = 4.64872629 # Bohr
params['tmp_dir']     = '../tmp/' # make a dir
# ph
# params['high_path']   = np.array([[[0.0000,0.0000,0.0000], [0.0000,0.5000,0.0000]],
# 	                              [[0.0000,0.5000,0.0000], [0.3333,0.3333,0.0000]],
# 								  [[0.3333,0.3333,0.0000], [0.0000,0.0000,0.0000]]])
# ele
params['high_path']   = np.array([[[0.0000,0.5000,0.0000], [0.3333,0.3333,0.0000]],
	                              [[0.3333,0.3333,0.0000], [0.0000,0.0000,0.0000]],
								  [[0.0000,0.0000,0.0000], [-0.3333,-0.3333,0.0000]],
								  [[-0.3333,-0.3333,0.0000], [0.0000,-0.5000,0.0000]]])
# mode = ph_sheng and ph_dist
# get these files from ShengBTE with a harmonic
params['omega_sheng'] = '../BTE.omega'
params['qp_full']     = '../BTE.qpoints_full'
params['qp']          = '../BTE.qpoints'

# mode = ph_dist, ele_dist, and ele_band
params['cdyna_fname'] = '../run_4_300K_keep/graphene_cdyna.h5'
params['snap_number'] = 1 # time, start 0
params['scale']       = 1000

# mode = ele_dist, ele_band
params['tet_fname']   = '../graphene_tet.h5' 
params['band_index']  = 0 # start 0

#-------------------------------------------------------------------------------------
'''
	read phonon distributions at snap_number and band_index
	return dist_func
'''
def ReadDist(params):

	snap_number = params['snap_number']
	band_index  = params['band_index']
	cdyna       = h5py.File(params['cdyna_fname'], 'r') # get the data 
	if params['mode'] == 'ph_dist':
		dist_func   = np.array(cdyna['ph_dyn_run_1']['snap_t_'+str(snap_number)])
		dist_func_0 = np.array(cdyna['ph_dyn_run_1']['snap_t_'+str(0)])
	else:
		dist_func   = np.array(cdyna['dynamics_run_1']['snap_t_'+str(snap_number)][:,band_index])
		dist_func_0 = np.array(cdyna['dynamics_run_1']['snap_t_'+str(0)][:,band_index])

	cdyna.close()

	if snap_number == 0:
		return dist_func_0
	else:
		return dist_func - dist_func_0

#-------------------------------------------------------------------------------------
'''
	distance
	time-consuming part
'''
def point2fixedAxis(point, fixedAxis):
    vector1 = point - fixedAxis[0]
    vector2 = point - fixedAxis[1]
    vector3 = fixedAxis[1] - fixedAxis[0]

    #k  = np.dot(fixedAxis[0] - point, fixedAxis[1] - fixedAxis[0])
    #k /= -np.square(np.linalg.norm(vector3))
    #dropFoot = k * (vector3) + fixedAxis[0]

    d = np.linalg.norm(np.cross(vector1, vector2)) / np.linalg.norm(vector3)

    return d

#-------------------------------------------------------------------------------------
'''
	choose the kpts at the high path
	return iqpt_index, iqpt_cho, RecVecs, omega
'''
def ChoosePts(params):

	high_path = params['high_path']
	band_ind  = params['band_index']
	alat      = params['alat']
	at        = params['at']
	Bohr2Ang  = 5.29177208590000E-01
	# RecVecs = 2*np.pi*RecVecs
	RecVecs   = np.linalg.inv((Bohr2Ang*alat*at).T)

	if params['mode'] in ['ph_sheng', 'ph_dist']:
		all_pts = np.loadtxt(params['qp_full'])
		omega   = np.loadtxt(params['omega_sheng']) * 0.6582122269864754
		num_pts = all_pts.shape[0]
	else:
		cdyna     = h5py.File(params['cdyna_fname'], 'r') # get the data 
		ryd2ev    = cdyna['band_structure_ryd'].attrs['ryd2ev'] 
		energy_ev = cdyna['band_structure_ryd'][:,band_ind] * ryd2ev  # eV
		cdyna.close()
		tet       = h5py.File(params['tet_fname'], 'r') 
		all_pts   = np.array(tet['kpts_all_crys_coord'])
		num_pts   = np.array(tet['num_kpts'])		
		tet.close()

	ipt_ind  = [] # start 0
	ipt_cho  = []
	ieig     = []

	# read qpt, kpt and index
	read = True
	qORk = 'q' if 'ph' in params['mode'] else 'e'
	for j in range(high_path.shape[0]):
		if os.path.isfile(params['tmp_dir'] + '%s%s_ind.npy' %(qORk, j)):
			ipt_ind.append(np.load(params['tmp_dir'] + '/%s%s_ind.npy' %(qORk, j)))
			read = False
		if os.path.isfile(params['tmp_dir'] + '%s%s_ipt.npy' %(qORk, j)):
			ipt_cho.append(np.load(params['tmp_dir'] + '%s%s_ipt.npy' %(qORk, j)))
			read = False
		
		if read:
			tmp_ind  = []; tmp_ipt = []
			path_dis = np.linalg.norm(high_path[j][1]-high_path[j][0])

			for i in range(num_pts):
				ipt    = all_pts[i, 2:] if 'ph' in params['mode'] else all_pts[i]
				ipt_p1 = np.linalg.norm(ipt-high_path[j][0])
				ipt_p2 = np.linalg.norm(ipt-high_path[j][1])
				dis  = point2fixedAxis(point=ipt, fixedAxis=high_path[j])
			
				if (dis < 0.001) and (ipt_p2 <= path_dis) and (ipt_p1 <= path_dis):
					tmp_ind.append(i)
					if 'ph' in params['mode']:
						tmp_ipt.append(all_pts[i, 2:])
					else:
						tmp_ipt.append(all_pts[i])
			
			tmp_ind = np.array(tmp_ind, dtype=int)
			tmp_ipt = np.array(tmp_ipt)

			# save 
			np.save(params['tmp_dir'] + "%s%s_ind.npy" %(qORk, j), tmp_ind)
			np.save(params['tmp_dir'] + "%s%s_ipt.npy" %(qORk, j), tmp_ipt)

			ipt_ind.append(tmp_ind)
			ipt_cho.append(tmp_ipt)

	# 
	for j in range(high_path.shape[0]):
		tmp_eig = []
		
		if 'ph' in params['mode']:
			for i in range(ipt_ind[j].shape[0]):
				tmp_eig.append(omega[int(all_pts[ipt_ind[j][i], 1]-1)])
		else:
			tmp_eig.append(energy_ev[ ipt_ind[j] ])

		tmp_eig = np.array(tmp_eig)
		ieig.append(tmp_eig)

	# print
	for i in range(high_path.shape[0]):
		print('    len(path-%s) = ' %(i+1), ipt_cho[i].shape)
	
	return ipt_ind, ipt_cho, RecVecs, ieig

#-------------------------------------------------------------------------------------
'''
	get x-axis in band
	return dis (x-axis in band)
'''
def GetX(ipt, high_path, RecVecs):

	for i in range(high_path.shape[0]):
		ipt[i] = np.dot(ipt[i], RecVecs)

	dis = []; dis_0_1 = 0
	for path in range(high_path.shape[0]):
		path_0 = np.dot(high_path[path][0], RecVecs)
		path_1 = np.dot(high_path[path][1], RecVecs)
		tmp    = []
		for p in ipt[path]:
			dis_ = np.linalg.norm(p - path_0)
			x    = dis_0_1 + dis_
			tmp.append(x)
		dis.append(np.array(tmp))
		dis_0_1 += np.linalg.norm(path_1 - path_0)

	return dis

#-------------------------------------------------------------------------------------
'''
	plot band and distributions
'''
def PlotDist(params, save_=False):

	high_path   = params['high_path']
	snap_number = params['snap_number']
	scale       = params['scale']
	band_index  = params['band_index']

	# get distributions 
	dist_func = ReadDist(params)

	# get ikpt_index, ikpt_cho, RecVecs
	ipt_index, ipt_cho, RecVecs, ieig = ChoosePts(params)

	# get dis
	dis = GetX(ipt_cho, high_path, RecVecs)

	# 
	xlim = [0]; x_0 = 0
	for path in range(high_path.shape[0]):
		path_0 = np.dot(high_path[path][0], RecVecs)
		path_1 = np.dot(high_path[path][1], RecVecs)
		x_0   += np.linalg.norm(path_1 - path_0)
		xlim.append(x_0)

	# plot
	plt.rc('font',family='Times New Roman')
	mpl.rcParams['xtick.direction'] = 'in'
	mpl.rcParams['ytick.direction'] = 'in'

	plt.figure(figsize=(6,5))
	plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)

	for i in range(high_path.shape[0]):
		if params["mode"] == 'ph_dist':
			ie   = ieig[i]
			dist = dist_func[ipt_index[i]]
			dist[np.where(dist_func[ipt_index[i]]<0)] = 0.0
			plt.scatter(np.repeat(dis[i], ie.shape[1]), ie, s=dist*scale)
			
		else:
			ie   = ieig[i][0]
			dist = dist_func[ipt_index[i]]
			plt.scatter(dis[i], ie, s=dist*scale)
		print(np.where(dist_func[ipt_index[i]]<0))
		plt.plot(dis[i], ie, color='grey') 
		

	if params['mode'] == 'ph_dist':
		plt.xlabel('q', size=20) 
		plt.ylabel('Frequency (meV)', size=20)
	else:
		plt.xlabel('k', size=20) 
		plt.ylabel('Energy (eV)', size=20)		

	#plt.ylim([-4.63, -3.5])
	plt.xlim([0, xlim[-1]])
	plt.show()

	if save_:
		if params['mode'] == 'ph_dist':
			plt.savefig('ph-%s.png' %snap_number, dpi=600)
		else:
			plt.savefig('ele-%s-%s.png' %(snap_number, band_index), dpi=600)

#-------------------------------------------------------------------------------------
def PlotPhBands(params, save_=False):

	high_path = params['high_path']

	# get ipt_index, ipt_cho, RecVecs
	ipt_ind, ipt_cho, RecVecs, ieig = ChoosePts(params)

	# get dis
	dis = GetX(ipt_cho, high_path, RecVecs)	

	xlim = [0]; x_0 = 0
	for path in range(high_path.shape[0]):
		path_0 = np.dot(high_path[path][0], RecVecs)
		path_1 = np.dot(high_path[path][1], RecVecs)
		x_0   += np.linalg.norm(path_1 - path_0)
		xlim.append(x_0)

	# plot
	plt.rc('font',family='Times New Roman')
	mpl.rcParams['xtick.direction'] = 'in'
	mpl.rcParams['ytick.direction'] = 'in'

	plt.figure(figsize=(6,5))
	plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)

	for i in range(high_path.shape[0]):
		if 'ele' in params['mode']:
			ie = ieig[i][0]
		else:
			ie = ieig[i]
		plt.plot(dis[i], ie, color='black') 

	if params['mode'] == 'ph_sheng':
		plt.xlabel('q', size=20) 
		plt.ylabel('Frequency (meV)', size=20) 
	else:
		plt.xlabel('k', size=20)
		plt.ylabel('Energy (eV)', size=20)

	#plt.ylim([-4.63, -3.5])
	plt.xlim([0, xlim[-1]])
	
	plt.show()

	if save_:
		if params['mode'] == 'ph_sheng':
			plt.savefig('ph_sheng.png', dpi=600)
			# save data

		else:
			plt.savefig('ele_band.png', dpi=600)
			# save data

#-------------------------------------------------------------------------------------
def main(params):

	print('========================< Plot >========================')

	# file check
	# mode = ph_sheng, ph_dist, and ele_dist
	if params['mode'] in ['ph_sheng', 'ph_dist']:
		assert os.path.isfile(params['omega_sheng']), exit('No %s! Exitting...' \
			                                          %params['omega_sheng'])
		assert os.path.isfile(params['qp']), exit('No %s! Exitting...' \
			                                          %params['qp'])
		assert os.path.isfile(params['qp_full']), exit('No %s! Exitting...' \
			                                          %params['qp_full'])
	elif params['mode'] in ['ph_dist', 'ele_dist', 'ele_band']:
		assert os.path.isfile(params['cdyna_fname']), exit('No %s! Exitting...'\
			                                          %params['cdyna_fname'])
	elif params['mode'] in ['ele_dist', 'ele_band']:
		assert os.path.isfile(params['tet_fname']), exit('No %s! Exitting...' \
			                                          %params['tet_fname'])
	else:
		exit('Unknow mode! mode = %s ??? Exiting...' %params['mode'])

	# plot ph_sheng
	if params['mode'] in ['ph_sheng', 'ele_band']:
		print('--> 1# info: plotting phonon dispersion with ShengBTE outs...')
		PlotPhBands(params)
	else:
		print('--> 1# info: plotting distribution at the snap time...')
		PlotDist(params)
	
	print('========================< Done >========================')

main(params)
