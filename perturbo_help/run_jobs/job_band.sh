#!/bin/bash
#SBATCH -p D2Part_test
#SBATCH -N 1
#SBATCH -n 56

# intel2020
source /es01/software/profile.d/intel2020.sh

# hdf5 1.12.0
export PATH=$PATH:/es01/software/apps/hdf5-1.12.0

# QE7.0
export PATH=$PATH:/es01/yeesuan/yeesuan623/soft/qe_7.0/bin

# wannier90
#export PATH=$PATH:/es01/yeesuan/yeesuan623/soft/qe_7.0/wannier90-3.0.0

# perturbo2.0_pristine
#export PATH=$PATH:/es01/yeesuan/yeesuan623/soft/qe_7.0/perturbo_pri

# perturbo_modification
#export PATH=$PATH:/es01/yeesuan/yeesuan623/soft/qe_7.0/perturbo_modi

#mpirun pw.x < scf.in > scf.log

mpirun pw.x < 2_cal_band.in > 2_cal_band.log &&

mpirun bands.x -i 3_plot_band.in > 3_plot_band.log &&

mpirun projwfc.x -i 4_plot_pband.in > 4_plot_pband.log &&

mpirun pp.x < 5_ppwf.in > 5_ppwf.log


