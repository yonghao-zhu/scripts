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
import linecache
import os, math
import matplotlib as mpl
import matplotlib.pyplot as plt
####################################
'''
    input file: XDATCAR
    output    : Slide.dat
'''
param = {}
param['TorF']     = True # True or False
param['mode']     = 'z' # xy1/2-->slide or z-->interlayer distance
# xy1 mode
param['the_atom'] = 20 # W11, start 1,
# xy2 mode 
param['the_atoms']= [2, 36]
param 
# z mode
param['layers']   = 4
param["layers_atoms"] = [32, 32, 12, 12] # from bottom to up
param['dis_two_layers'] = [1,2] # start 1
param['plot_']    = True
param['y_range']  = [2, 4, 0.5]
param['x_range']  = [0, 2000, 500]
####################################
def ReadXdatcar(filename, the_atom=0):

    if os.path.isfile('xdatcar.npy'):
        xdatcar_C = np.load('xdatcar.npy')

        the_atom_pos = xdatcar_C[:,the_atom-1,:]

    else:
        lattice = np.mat([[float(linecache.getline(filename,i+1).split()[0]),
            float(linecache.getline(filename,i+1).split()[1]),
            float(linecache.getline(filename,i+1).split()[2])] for i in range(2,5)])

        seventh_line = linecache.getline(filename,7).split()
        
        atoms = int(sum([float(i) for i in seventh_line]))
        
        comments = ['Direct configuration=']
        
        for i in range(7):
            comments.append(linecache.getline(filename,i+1))
        
        xdatcar_D = np.loadtxt(filename,comments=comments).reshape(-1,atoms,3)
        xdatcar_C =np.array([np.dot(i,lattice) for i in xdatcar_D])

        the_atom_pos = xdatcar_C[:,the_atom-1,:]

        np.save('xdatcar.npy', xdatcar_C)

    return the_atom_pos,xdatcar_C

#-----------------------------------

def plot(x, y, y_range, x_range):

    length_y = int((max(y_range[:2]) - min(y_range[:2])) / y_range[2])
    length_x = int((max(x_range[:2]) - min(x_range[:2])) / x_range[2])
    
    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'        

    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)

    plt.plot(x,y,c='black',zorder=0,lw=4)
    
    # y
    plt.ylabel("Distance (Ångström)",size=20)
    plt.ylim(y_range[:2])
    plt.yticks([round(y_range[0] + i*y_range[2], 2) for i in range(length_y+1)],
        [round(y_range[0] + i*y_range[2], 2) for i in range(length_y+1)],size=15)

    # x
    plt.xlabel("Time (fs)",size=20)
    plt.xlim(x_range[:2])
    plt.xticks([x_range[0] + i*x_range[2] for i in range(length_x+1)],
        [x_range[0] + i*x_range[2] for i in range(length_x+1)],size=15)        

    plt.show()

#-----------------------------------

def InterlayerSlide1(param):

    the_atom = param['the_atom']

    # read pos
    the_atom_pos,xdatcar_C = ReadXdatcar(filename='XDATCAR', the_atom=the_atom)
    print('the_atom_pos.shape =', the_atom_pos.shape)

    # calculate distance
    time_min_pos = the_atom_pos[0,:]
    distance = []
    for i in the_atom_pos:
        dis = math.sqrt( (i[0]-time_min_pos[0])**2 + 
                         (i[1]-time_min_pos[1])**2 )
        distance.append(dis)

    # write files
    filename = 'atom-%s-Slide.dat' %the_atom
    with open(filename, 'w+') as f:
        for time in range(len(distance)):
            f.writelines('  '+'%04d' % time+'  '+'%0.6f' %distance[time]+'\n')

    # plot
    if param['plot_']:
        plot(x=[i for i in range(len(distance))], y=distance, 
                y_range=param['y_range'], x_range=param['x_range'])

#-----------------------------------
def InterlayerSlide2(param):

    the_atoms = param["the_atoms"]

    # read pos1
    pos1,xdatcar_C = ReadXdatcar(filename='XDATCAR', the_atom=the_atoms[0])

    # read pos2 
    pos2,xdatcar_C = ReadXdatcar(filename='XDATCAR', the_atom=the_atoms[1])

    # calculate distance
    distance = []
    for i in range(pos1.shape[0]):
        dis = np.sqrt( (pos1[i, 0]-pos2[i, 0])**2 
                     + (pos1[i, 1]-pos2[i, 1])**2 )
        distance.append(dis)    

    # write files
    filename = 'atom-%s-%s-Slide.dat' %(the_atoms[0], the_atoms[1])
    with open(filename, 'w+') as f:
        for time in range(len(distance)):
            f.writelines('  '+'%04d' % time+'  '+'%0.6f' %distance[time]+'\n')

    # plot
    if param['plot_']:
        plot(x=[i for i in range(len(distance))], y=distance, 
                y_range=param['y_range'], x_range=param['x_range'])

#-----------------------------------

def InterlayerDistance(param):

    # read pos
    the_atom_pos,xdatcar_C = ReadXdatcar(filename='XDATCAR')
    print('xdatcar_C.shape =', xdatcar_C.shape)   

    layers_atoms = param["layers_atoms"]

    layers_index = []
    # create layers index
    for i in range(param['layers']):
    
        if i == 0:
            tmp = [0, layers_atoms[0]]

        else:
            start = sum(layers_atoms[:i])

            if i+1 == param["layers"]:
                end   = sum(layers_atoms[:])

            else:
                end   = sum(layers_atoms[:i+1])

            tmp = [start, end]

        layers_index.append(tmp)

    print("layer index=", layers_index)

    distance = []
    for time in range(xdatcar_C.shape[0]):

        tmp = sorted( list(xdatcar_C[time,:,2]) )  
        
        layer_1 = tmp[ layers_index[param['dis_two_layers'][0]-1][0] \
                      :layers_index[param['dis_two_layers'][0]-1][1] ]
        
        layer_2 = tmp[ layers_index[param['dis_two_layers'][1]-1][0] \
                      :layers_index[param['dis_two_layers'][1]-1][1] ]
        
        dis =   sum(layer_2)/layers_atoms[param['dis_two_layers'][1]-1] \
              - sum(layer_1)/layers_atoms[param['dis_two_layers'][0]-1]

        distance.append(dis)

    # write files
    filename = 'breathing.dat'
    with open(filename, 'w+') as f:
        for time in range(xdatcar_C.shape[0]):
            f.writelines('%04d' %time + '  ' + '%0.6f' %distance[time] + '\n')

    # plot
    if param['plot_']:
        print(distance[0])
        plot(x=[i for i in range(len(distance))], y=distance, 
                x_range=param['x_range'], y_range=param['y_range'])

####################################
if param['TorF']:
    if os.path.isfile('XDATCAR'):
        if param['mode'] == 'xy1':
            print('xy mode Running.....')
            InterlayerSlide1(param)

        if param["mode"] == 'xy2':
            print('xy mode Running.....')
            InterlayerSlide2(param)

        if param['mode'] == 'z':
            print('z mode Running......')
            InterlayerDistance(param)
    else:
        print('No XDATCAR! Exiting...')

