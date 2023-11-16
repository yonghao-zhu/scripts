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
    plot the etot form MOVEMENT
'''

import numpy as np
import os

import matplotlib as mpl
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------
def GetDataFromMment(movement, type_, research):

	if not os.path.isfile(movement):
		exit('No MOVEMENT! Exiting...')

	with open(movement, 'r') as f:
		data_row = []
		data_typ = []
		for il in f.readlines():
			if research in il:
				data_row.append(il)
				if type_ == 'Etot':
					line = il.split()
					data_typ.append( float(line[6]) )
	return data_typ

#-------------------------------------------------------------------------
def PlotData(data):

	plt.rc('font',family='Times New Roman')
	mpl.rcParams['xtick.direction'] = 'in'
	mpl.rcParams['ytick.direction'] = 'in'

	plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.18)

	steps = [i for i in range(len(data))]
	data  = [i-data[0] for i in data]

	plt.plot(steps[-40:], data[-40:])

	plt.show()

#-------------------------------------------------------------------------
def main():

	movement = 'MOVEMENT'
	type_    = 'Etot'
	research = 'HSE_Iter'
	data     = GetDataFromMment(movement, type_, research)
	PlotData(data)

if __name__ == '__main__':
	main()