#!/bin/bash
#SBATCH -p vip_07
#SBATCH -N 2
#SBATCH -n 24

source /public3/soft/modules/module.sh

module load intel/19.0.3-bscc-public3

# QE6.5, including wannier90.x and epw.x
#module load QE/6.5-intel19-public3

# QE7.0
export PATH=$PATH:/public3/home/sc30334/works/zhuyh/soft/QE_7.0/bin

# test run
#srun pw.x < DISP.BAs_sc.in.001 > out_1.log

end_ph=140
run_start=79 # generally fixed
prefix="DISP.BAs_sc.in."

#run jobs
#-------------------------------------------------------------
echo Job starts at `date`
#-------------------------------------------------------------

for i in `seq $run_start $end_ph`
do
	printf "\n=======================================\n"
	printf "Running at $i...... \n"
	
	# job in.xxx
	scf_path="./force_run/$prefix$i"  &&

	# check ph_x
	if [ ! -d ./$scf_path ]; then exit 0; fi	
	if [ -f ./$scf_path/out.log ]; then exit 0; fi

	# run job
	cd $scf_path &&
	echo "run $scf_path..." &&
	
	srun pw.x < scf.in > out.log &&
	
	if [ -d ./tmp ]; then rm -r tmp; fi
	
	cd ../../ &&
	
	sleep 1
done

#-------------------------------------------------------------
echo Job ends at `date`
#-------------------------------------------------------------
