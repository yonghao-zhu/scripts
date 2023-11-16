#!/bin/sh
#SBATCH --exclusive   
#SBATCH --partition=normal
#SBATCH --job-name=GaN
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --gres=gpu:4  
#SBATCH --exclude=gn25,gn47

module load intel/2016
module load cuda/8.0
source /data/home/zhuyh/codes/PWMAT/20230726/pwmat.env

#mpirun -np $SLURM_NPROCS -rdma PWmat -host 10.3.0.100 50001 | tee output

#############################################################
START_T=0  # fs
STEP_T=10  # fs
END_T=1000 # fs

for time in `seq $START_T $STEP_T $END_T`
do
	SUFX=$( printf "%04d" "$time" )
	cd ./RunDos/"dos-${SUFX}"
	
	# total dos
	if [ ! -d ./total_dos ]; then mkdir total_dos; fi
	if [ ! -f ./total_dos/bpsiiofil100001 ]; then
		cd total_dos &&
		# copy files
		if [ ! -f ./N.SG15.PBE.UPF ]; then ln -sf ../../../N.SG15.PBE.UPF; fi
		if [ ! -f ./Ga.SG15.PBE.UPF ]; then ln -sf ../../../Ga.SG15.PBE.UPF; fi
		if [ ! -f ./atom.config ]; then ln -sf ../atom.config; fi
		if [ ! -f ./IN.WG ]; then ln -sf ../IN.WG; fi
		if [ ! -f ./OUT.EIGEN ]; then ln -sf ../OUT.EIGEN; fi
		if [ ! -f ./etot.input ]; then cp ../../../etot-dos.input etot.input; fi
		mpirun -np $SLURM_NPROCS -rdma PWmat -host 10.3.0.100 50001 | tee output &&
		cd ../
	fi
	sleep 1
	
	# occ dos
	if [ ! -d ./occ_dos ]; then mkdir occ_dos; fi
	if [ ! -f ./occ_dos/bpsiiofil100001 ]; then
		cd occ_dos &&
		# copy files
		if [ ! -f ./N.SG15.PBE.UPF ]; then ln -sf ../../../N.SG15.PBE.UPF; fi
		if [ ! -f ./Ga.SG15.PBE.UPF ]; then ln -sf ../../../Ga.SG15.PBE.UPF; fi
		if [ ! -f ./atom.config ]; then ln -sf ../atom.config; fi
		if [ ! -f ./IN.WG ]; then ln -sf ../IN.WG; fi
		if [ ! -f ./OUT.EIGEN ]; then ln -sf ../OUT.EIGEN; fi
		if [ ! -f ./IN.OCC_ADIA ]; then ln -sf ../IN.OCC_ADIA; fi
		if [ ! -f ./etot.input ]; then cp ../../../etot-dos.input etot.input; fi
		echo "IN.OCC_ADIA=T" >> etot.input
		mpirun -np $SLURM_NPROCS -rdma PWmat -host 10.3.0.100 50001 | tee output &&
		cd ../
	fi
	sleep 1

	cd ../../
	
done
