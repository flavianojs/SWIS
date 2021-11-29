#!/bin/bash

	groundstate=0
	    compile=1
	  executing=1
forceexecuting=0
	    lattice=0
	 occupation=0
dispersionplot=0

export OMP_NUM_THREADS=4
# ulimit -s unlimited

source /usr/local/bin/compilervars.sh intel64
# source /usr/local/bin/compilervars.sh intel64
# source /usr/local/intel/mkl/bin/mklvars.sh intel64

XX=$1
shift
# Get optional flags
while (( "$#" )); do
   case $1 in
      (uff|iff|osx|juropa|juropatest|jureca|juqueen)
        addplatform="$1"
        platform="PLATFORM=$1"
        shift
        ;;
      (debug)
        adddebug=${adddebug}"_$1"
        debug="DEBUG=$1"
        shift
        ;;
      (cleandep|cleandebug|cleanobj|cleanmod|cleanexe|clean|cleanall|recompile)
        rule="$1 ${rule}"
        shift
        ;;
      (*.exe)
        if [[ -z "${1:0:${#1}-4}" ]] ; then
          echo "Empty filename!"
          exit 1
        fi
        filename=$(echo FILE=$1)
        shift
        ;;
      (verbose)
        verbose="--debug=$1"
        shift
        ;;
      (run)
        forceexecuting=1
        shift
        ;;
      (*)
        echo "Illegal option: $1"
        exit 1
   esac
done


#--- Ground state calculation ------------------------------------------

	if [ $groundstate -eq 1 ]; then
		cp inputfiles/lattice_$XX\.datbkp inputfiles/lattice_$XX\.dat 
		cp inputfiles/lattice_$XX\.dat ground\ state\ program/lattice_sky.dat 
		cp inputcard_$XX\.inp ground\ state\ program/inputcard_sky.inp 
		cd ground\ state\ program/

		echo "Ground state program starting:"
		./topo_sw.x

		cp plot.dat ../inputfiles/lattice_$XX\.dat 
		cd .. 
	fi

#--- Dispersion calculation --------------------------------------------
	#--- Compiling ------------------------------------------------------

	if [ $compile -eq 1 ]; then 
		cd ~/Dropbox/scripts/dispersion-program
		make $rule $platform $debug $filename $verbose
		cd -
	fi

	if [[ $rule =~ "clean" ]]; then 
		exit
	fi

	#--- Executing ------------------------------------------------------

	# source /usr/local/bin/compilervars.sh intel64

	runsuccess=true
	if [ $executing -eq 1 ]; then

		# Did the executable change
      if ! cmp ~/Dropbox/scripts/dispersion-program/bin/main.exe ~/Dropbox/scripts/dispersion-program/bin/.main.exe >/dev/null 2>&1
      then
         forceexecuting=1
      fi

		if [ $forceexecuting -eq 1 ]
		then
			echo -e "\nSpinwave dispersion program starting because of recompilation:"
			~/Dropbox/scripts/dispersion-program/bin/main.exe  inputcard_$XX\.inp
			# ~/Dropbox/scripts/dispersion-program/bin/main_y.exe  inputcard_$XX\.inp
         if [ "$?" -eq "0" ]; then
				mv precession_$XX.dat disp_analy_$XX.dat disp_unfol_$XX.dat  dispersion_$XX.dat latticExt_$XX.dat kpath_$XX.dat outputfiles/

   			cp ~/Dropbox/scripts/dispersion-program/bin/main.exe ~/Dropbox/scripts/dispersion-program/bin/.main.exe
   			cp inputcard_$XX\.inp .inputcard_$XX\.inp
			else
			  runsuccess=false
			fi
		else
			 # Did the inputcard change?
			if ! cmp inputcard_$XX\.inp .inputcard_$XX\.inp >/dev/null 2>&1
			then
				echo -e "\nSpinwave dispersion program starting:"
				~/Dropbox/scripts/dispersion-program/bin/main.exe  inputcard_$XX\.inp
			# ~/Dropbox/scripts/dispersion-program/bin/main_y.exe  inputcard_$XX\.inp
				if [ "$?" -eq "0" ]; then
					mv precession_$XX.dat disp_analy_$XX.dat disp_unfol_$XX.dat  dispersion_$XX.dat latticExt_$XX.dat kpath_$XX.dat outputfiles/
   				cp inputcard_$XX\.inp .inputcard_$XX\.inp
	  			else
					runsuccess=false
				fi
			else
				echo "Calculated data are up to date! Plotting only..."
			fi # Did the inputcard change?
		fi # Did the executable change

	fi	# executing == 1

#--- Plottings ---------------------------------------------------------
	#--- Lattice Plotting -----------------------------------------------

	if [ $lattice -eq 1 ]; then
		# basisname=$(grep "basisname" inputcard_$XX.inp | awk '{print $3}')
		# echo $basisname
		# gr lattice_sky.py $basisname

		# gr lattice_sky.py inputfiles/spinconfig_02.txt
		# gr lattice_sky.py "inputfiles/spinconfig_AF20x1.txt" "inputfiles/lattice_AF20x1.dat"
		# gr lattice_sky.py inputfiles/spinconfig_AF20x20AS_spiral.txt inputfiles/lattice_$XX.dat

		gr lattice_sky.py inputfiles/spinconfig_$XX.txt inputfiles/lattice_$XX.dat spinconfig_$XX.png
		if [ "$?" -eq "0" ]; then open spinconfig_$XX.png; fi
		# open lattice_sky.html
	fi

	#--- Occupation number ----------------------------------------------

	if [ $occupation -eq 1 ]; then
		gnuplot occupation.gnu
		open occupation.html
	fi	

	#--- Dispersion -----------------------------------------------------
	if [ $dispersionplot -eq 1 ] && $runsuccess ; then
		gnuscript="disp_unfol_multiPNG.gnu"

		maxomega=$(grep "maxomega" inputcard_$XX.inp | awk '{print $3}')
		minomega=$(grep "minomega" inputcard_$XX.inp | awk '{print $7}')

		gnuplot -e "specifier='$XX'; ymax='$maxomega'; ymin='$minomega'" $gnuscript
		if [ "$?" -eq "0" ]; then open disp_unfol_$XX.png ; fi
		# open disp_unfol.png 
			
	fi
