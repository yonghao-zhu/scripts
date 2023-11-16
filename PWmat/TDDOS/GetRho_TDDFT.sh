#! /bin/bash

source ~/.bashrc

# 1: gen RHO; 2: gen xsf; 3: plot dos
run_mode=2 # 1, 2, 3

START_T=0  # fs
STEP_T=40  # fs
END_T=1200 # fs

for time in `seq $START_T $STEP_T $END_T`
do
	SUFX=$( printf "%04d" "$time" )
	cd ./RunDos/"dos-${SUFX}"
	
	sleep 1
	
	if [ $run_mode -eq 1 ]; then
		echo `pwd`
		if [ ! -d ./charge_density ]; then mkdir charge_density; fi
		cd charge_density
		sed "s/X/$time/g" ../../../kband.txt > ./kband.txt

		ln -sf ../atom.config
		ln -sf ../IN.WG OUT.WG
		ln -sf ../../../OUT.TDDFT1
		ln -sf ../TDDFT.RHO OUT.RHO
		ln -sf ../../../OUT.GKK
		ln -sf ../../../REPORT

		~/scripts/split_tddft_occ_epfix.x &&
		if [ ! -d TDDOS ]; then mkdir TDDOS; fi
		cp atom.config TDDOS/restart.config &&
		~/scripts/generate_defect_rho.x &&

		rm atom.config OUT.WG OUT.TDDFT1 OUT.RHO OUT.GKK REPORT &&
		cd ..
		
	elif [ $run_mode -eq 2 ]; then
		echo `pwd`
		cd charge_density
		file_count=$(ls -l | grep 'RHO1' | wc -l)
		for ((i=1; i<=$file_count; i++)); do
			rho_sufx=1-$i 
			if [ ! -f ./REPORT ]; then ln -sf ../../../REPORT; fi
			if [ ! -f ./atom.config ]; then ln -sf ../atom.config; fi
			convert_rho.x IN.RHO$rho_sufx &&
			mv RHO.xsf RHO.$rho_sufx.xsf &&
			rm REPORT atom.config 
		done
		cd ..

	elif [ $run_mode -eq 3 ]; then
		echo `pwd` &&
	
		if [ -f ./total_dos/bpsiiofil100001 ]; then
			cd ./total_dos/
			#ln -sf ../../../OUT.FERMI &&
			#sed "s/X/0/g" ../../../DOS.input > ./DOS.input &&
			#plot_DOS.py &&
			#plot_DOS_interp.x &&
			#plot_DOS.py &&

			rm *UPF OUT.EIGEN IN.WG IN.OCC_ADIA OUT.FERMI OUT.IND_EXT_KPT OUT.KPT OUT.NONSCF OUT.TDDFT1 
			rm TIMELOG TDDFT.RHO OUT.GKK OUT.SYMM OUT.WG ORIGIN.INDEX OUT.RHO
			cd ../
		fi 

		if [ -f ./occ_dos/bpsiiofil100001 ]; then
			cd ./occ_dos/
			#ln -sf ../../../OUT.FERMI &&
			#sed "s/X/1/g" ../../../DOS.input > ./DOS.input &&
			#plot_DOS.py &&
			#plot_DOS_interp.x &&
			#plot_DOS.py &&

			rm *UPF OUT.EIGEN IN.WG IN.OCC_ADIA OUT.FERMI OUT.IND_EXT_KPT OUT.KPT OUT.NONSCF OUT.TDDFT1 
			rm TIMELOG TDDFT.RHO OUT.GKK OUT.SYMM OUT.WG ORIGIN.INDEX OUT.RHO
			cd ../
		fi
		
	fi # end run_mode
	
	cd ../../
	
done

