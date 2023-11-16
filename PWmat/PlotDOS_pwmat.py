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
	plot dos with pwmat format
	input: DOS.totalspin or DOS.totalspin_projected (DOS)
		   OUT.OCC (SCF)
	NOTE: VBMs set to zero!
	author: yonghao_zhu@163.com
''' 

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt 
import os, linecache 

#-------------------------------------------------------------------------------
def ReadOCC(occ_fname):
	# get nkpts and nbands
	nkpts = 0; KPOINTS_line = []
	with open(occ_fname, 'r') as f:
		lines  = f.readlines()
		nlines = len(lines)
		for i in range(nlines):
			if 'KPOINTS' in lines[i]:
				nkpts += 1
				KPOINTS_line.append(i)
	print('nkpts  = ', nkpts)
	if nkpts != 1:
		nbands = KPOINTS_line[1] - KPOINTS_line[0] - 2
	else:
		nbands = nlines - 2
	print('nbands = ', nbands)
	# get occ
	occ = np.zeros([nkpts, nbands, 3])
	for ik in range(nkpts):
		for ib in range(nbands):
			iline = linecache.getline(occ_fname, ik*(nbands+2)+ib+3).split()
			occ[ik, ib, 0] = int(iline[0])
			occ[ik, ib, 1] = float(iline[1])
			occ[ik, ib, 2] = float(iline[2])
	return occ

def ReadDOS(doscar):
	# read Line 1
	prefix = linecache.getline(doscar, 1).split()
	print(prefix)
	# read dos
	dos = np.loadtxt(doscar)
	print('dos.shape = ', dos.shape)

	return dos, prefix

def PlotDOS(occ, dos, prefix, ene_rang, dos_rang):
	# get VBM
	vbm = np.max(occ[np.where(occ[:,:,2] > 0.0)][:, 1])
	print('vbm = ', vbm)
	# plot
	plt.rc('font',family='Times New Roman')
	mpl.rcParams['xtick.direction'] = 'in'
	mpl.rcParams['ytick.direction'] = 'in'

	plt.figure(figsize=(6,5))
	plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)

	for i in range(dos.shape[1] - 2):
		plt.plot(dos[:, 0]-vbm, dos[:, i+2], lw=4) 
	plt.legend(prefix[3:], fontsize='large')
	plt.xlabel('Energy (eV)', size=20) 
	plt.ylabel(f'DOS (eV⁻¹)', size=20) 
	
	plt.xticks(size=20)
	plt.yticks(size=20)
	# plot VBM
	plt.axvline(x=0.0,lw=2)
	plt.xlim(ene_rang)
	plt.ylim(dos_rang)


	plt.show()

#------------------------[1:]-------------------------------------------------------
def main():
	params = {}
	params['doscar']  = './DOS.totalspin'
	params['out.occ'] = '../2_scf/OUT.OCC'
	params['save_']   = False # True or False
	occ        = ReadOCC(params['out.occ'])
	dos,prefix = ReadDOS(params['doscar'])
	ene_rang   = [-1.5, 2] 
	dos_rang   = [0, 3]
	PlotDOS(occ, dos, prefix, ene_rang, dos_rang)

if __name__ == '__main__':
	main()