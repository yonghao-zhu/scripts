#!/bin/bash
#SBATCH -p D2Part_test
#SBATCH -N 1
#SBATCH -n 36

# intel2020
source /es01/software/profile.d/intel2020.sh

# hdf5 1.12.0
#export PATH=$PATH:/es01/software/apps/hdf5-1.12.0

# QE7.0
export PATH=$PATH:/es01/yeesuan/yeesuan623/soft/qe_7.0/bin

# wannier90
#export PATH=$PATH:/es01/yeesuan/yeesuan623/soft/qe_7.0/wannier90-3.0.0

# perturbo2.0_pristine
#export PATH=$PATH:/es01/yeesuan/yeesuan623/soft/qe_7.0/perturbo_pri

# perturbo_modification
#export PATH=$PATH:/es01/yeesuan/yeesuan623/soft/qe_7.0/perturbo_modi

#mpirun pw.x < scf.in > scf.log

#mpirun pw.x < nscf.in > nscf.log

start_ph=1 # fixed
end_ph=29
run_start=2 # generally fixed
prefix="bas"
CASE=ph

#run jobs
#-------------------------------------------------------------
echo Job starts at `date`
#-------------------------------------------------------------

# mkdir save
if [ ! -d ./save/"$prefix.phsave" ]; then mkdir -p ./save/"$prefix.phsave"; fi &&

for i in `seq $run_start $end_ph`
do
	printf "\n=======================================\n"
	printf "Running at ph-$i...... \n"
	
	# ph_x name
	ph_x="ph-$i"  &&

	# check ph_x
	if [ ! -d ./$ph_x ]; then exit 0; fi
	
	# run in ph_x
	cd ./$ph_x &&
	mkdir ./tmp &&
	cd ./tmp/ &&
	ln -sf ../../tmp/"$prefix.save" &&
	cd .. &&
	mpirun ph.x $PH < ${CASE}.in > ${CASE}.out &&

	sleep 1 &&
	# remove wfc files
	rm ./tmp/*wfc* &&

	dvscf_0="$prefix.dvscf1" &&
	dvscf_1="$prefix.dvscf_q$i" &&
	cp ./tmp/_ph0/"$prefix.phsave"/* ../save/"$prefix.phsave"/ &&
	cp ./*dyn* ../save

	cp ./tmp/_ph0/"$prefix.q_$i"/$dvscf_0 ../save/$dvscf_1 &&	
	# remove wfc files in prefix.q_x
	rm ./tmp/_ph0/"$prefix.q_$i"/*wfc* &&	
	cd .. &&

	sleep 1
done

#-------------------------------------------------------------
echo Job ends at `date`
#-------------------------------------------------------------


