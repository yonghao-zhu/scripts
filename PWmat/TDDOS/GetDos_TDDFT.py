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
    generate the TDDOS
    run me in the TDDFT dir
'''

import os
import linecache

#--------------------------------------------------------------------
def GenDosDir(fname_MOVEMENT, total_steps, md_step, out_step):

    # get RunDos dir
    if not os.path.isdir('RunDos'):
        print("mkdir RunDos...")
        os.system("mkdir RunDos")

    # get out_time
    md_time  = int(total_steps * md_step) # fs
    out_step = [i for i in range(0, md_time+out_step, out_step)]
    
    # read natoms from MOVEMENT
    if not os.path.isfile(fname_MOVEMENT):
        exit('No %s File! Exiting...' %fname_MOVEMENT)

    natoms = int(linecache.getline(fname_MOVEMENT, 1).split()[0])
    print('--> natoms = ', natoms)

    for it in out_step:
        # get dos_xxxx dir
        if not os.path.isdir("RunDos/dos-%04d" %it):
            print("mkdir RunDos/dos-%04d" %it)
            os.system("mkdir RunDos/dos-%04d" %it)

        # copy files
        digits = len(str(it))
        if it == 0:
            digits = 0
        suffix = "%0.6fE+0%d" %(it / (10**digits), digits)
        
        os.system("ln -sf ../../TDDOS/OUT.EIGEN.%s ./RunDos/dos-%04d/OUT.EIGEN" \
                    %(suffix, it))
        os.system("ln -sf ../../TDDOS/OUT.WG.%s ./RunDos/dos-%04d/IN.WG" \
                    %(suffix, it))
        os.system("ln -sf ../../TDDOS/OUT.OCC_ADIA.%s ./RunDos/dos-%04d/IN.OCC_ADIA" \
                    %(suffix, it))
        os.system("ln -sf ../../TDDOS/OUT.RHO.%s ./RunDos/dos-%04d/TDDFT.RHO" \
                    %(suffix, it))

        # get atom.config
        iIter        = int(it / md_step)
        print(iIter)
        it_ini_lines = iIter*(natoms*3 + 3 + 5 + 1) + 1

        with open('./RunDos/dos-%04d/atom.config' %it, 'w+') as f:
            for iatom in range(natoms*3+3+6):
                lines = linecache.getline(fname_MOVEMENT, it_ini_lines+iatom)
                f.writelines(lines)

#--------------------------------------------------------------------
def main():

    # Generate Dir to calculate tddos
    # total dirs = total_steps * md_step / out_step
    fname_MOVEMENT = 'MOVEMENT'
    total_steps    = 10000 
    md_step        = 0.1 # fs
    out_step       = 10  # fs
    GenDosDir(fname_MOVEMENT, total_steps, md_step, out_step)

if __name__ == "__main__":
    main()

