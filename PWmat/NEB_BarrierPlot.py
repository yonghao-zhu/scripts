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

import numpy as np
np.set_printoptions(precision=15)
import os
import linecache

import matplotlib as mpl
import matplotlib.pyplot as plt

#-------------------------------------------------------
def GetBarrier(barrier_file, imag):

	barrier = []
	with open(barrier_file, 'r') as f:
		iter_list = []
		lines = f.readlines()
		for i in range(len(lines)):
			if 'iter=' in lines[i]:
				iter_list.append(i)
		
		for i in range(imag):
			barrier.append(float(lines[iter_list[-1]+i+1].split()[1]))

	return barrier

def PlotBarrier(barrier_file, imag):

	barrier = np.array(GetBarrier(barrier_file, imag))

	plt.rc('font',family='Times New Roman')
	mpl.rcParams['xtick.direction'] = 'in'
	mpl.rcParams['ytick.direction'] = 'in'	
	plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.18)

	plt.xlabel('Step', size=25); plt.ylabel('Energy (eV)', size=25)
	plt.xlim(1, imag); plt.xticks(size=18)
	plt.ylim(-0.2, 1.2); plt.yticks(size=18)

	x = np.array([i for i in range(1, 1+imag)])
	plt.plot(x, barrier-barrier[0], 'o-', linewidth=4, markersize=16)

	print(barrier-barrier[0])
	print('energy barrier = ', max(barrier-barrier[0]), 'eV')

	plt.show()

#--------------------------------------------------------
def main():

	barrier_file = './NEB.BARRIER'
	imag         = 9 # total

	PlotBarrier(barrier_file, imag)

if __name__ == '__main__':
	main()

