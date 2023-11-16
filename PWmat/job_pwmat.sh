#!/bin/sh
#SBATCH --exclusive   
#SBATCH --partition=normal
#SBATCH --job-name=TEST
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gres=gpu:4  

module load intel/2016
module load cuda/8.0
source /data/home/zhuyh/codes/PWMAT/20230726/pwmat.env

mpirun -np $SLURM_NPROCS -rdma PWmat | tee output
