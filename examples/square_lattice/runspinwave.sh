#!/bin/bash
        compile=1
      executing=1
 forceexecuting=0
        lattice=0
     occupation=0
          kpath=1
 dispersionplot=1 
         spirit=1
     scale_pair=0

spirit_config_file='input_spirit.cfg'
spirit_initial_state='spirit_output/spinconfig_initial.ovf'

scale_input='inputfiles/pair.txt'
scale_output='inputfiles/pair_temp.txt'

gnuscript="disp_unfol_multi.gnu"

SWIS_path="/Users/flavianojs/Downloads/Spinwaves-dispersion"

#--- Setup environment ------------------------------------------

export OMP_NUM_THREADS=1
# ulimit -s unlimited

host=`hostname`
echo Hostname: $host

SWcode_executable="main_$host.exe"

if ! [ -f "$SWcode_executable" ]; then
    echo "$SWcode_executable does not exist. Forcing compilation."
    rule="recompile ${rule}"
fi

# For flaviano's MacBook Pro
if [ $host == 'Flavianos-MacBook-Pro.local' ] || [ $host == 'tsf-452-wpa-4-009.epfl.ch' ] ; then
    source /opt/intel/bin/compilervars.sh intel64
    source /opt/intel/mkl/bin/mklvars.sh intel64
    export OMP_NUM_THREADS=8

# For mb-santos
elif [ $host == 'mb-dossantos' ] ; then
    source /usr/local/bin/compilervars.sh intel64
    # source /usr/local/intel/mkl/bin/mklvars.sh intel64
    export OMP_NUM_THREADS=8

# For theospc47
elif [ $host == 'theospc47' ] ; then
	source /opt/intel/bin/compilervars.sh intel64
	source /opt/intel/mkl/bin/mklvars.sh intel64
    export OMP_NUM_THREADS=16
else
    echo "Unknown Hostcomputer: $host"
fi

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

    if [ $scale_pair -eq 1 ]; then
        echo -e '******* Scale pair ***********'
        python scale_pair.py $scale_input $scale_output
    fi

    if [ $spirit -eq 1 ]; then
        echo -e '******* Ground state ***********'

        # spirit_initial_state='spirit_output/spinconfig_Mn5Ge3_initial_random.ovf'
        spirit_final_state=$(grep "basisname" inputcard_$XX.inp | awk '{print $3}' | cut -d '"' -f 2 | cut -d ',' -f 2)

        should_execute=false
        # if [ $forceexecuting -eq 1 ]; then
        #     should_execute=true
        # fi

        pairfile=$(grep "pairfile" inputcard_$XX.inp | awk '{print $3}' | cut -d '"' -f 2 | cut -d ',' -f 2)

        declare -a files_to_check=(
                        $spirit_config_file
                        $pairfile
                        )
        # now loop through the above array
        for file in "${files_to_check[@]}"
        do
           if ! cmp --silent $file temp/$file
           then
               echo "Modified file:" $file
               should_execute=true
           fi
        done

        if $should_execute
        then
            echo -e "\nSpirit code starting:"
            # Testing if the Spirit inputcard has a simple unit cell of 1 1 1
            unit_cell_expansion=$(grep "n_basis_cells" $spirit_config_file | awk '{print $2 $3 $4 }' )

            if [ $unit_cell_expansion -ne "111" ]; then
                echo Set 'n_basis_cells' to 1 1 1 in the Spirit inputcard: $spirit_config_file
                exit
            fi

            python runspirit.py $spirit_config_file $spirit_initial_state $spirit_final_state

            #if the program has run successfully, make the copying and backups 
            if [ "$?" -eq "0" ]; then
                #Copying files to be check for modifications on the next run
                for file in "${files_to_check[@]}"
                do
                    echo $file
                    #Copy the file and the folder structuring
                    rsync -R $file temp/
                done
            fi
        else
            echo "Nothing to be done."
        fi
    fi

#--- Dispersion calculation --------------------------------------------
    #--- Compiling ------------------------------------------------------

    compiled=false
    if [ $compile -eq 1 ]; then 
        echo -e '\n******* Spin-wave code compilation ***********'
        cd $SWIS_path
        make $rule $platform $debug $filename $verbose | tee compilation_log.dat

        line=$(head -n 1 compilation_log.dat)
        rm compilation_log.dat
        cd - >/dev/null #not shows the path when going back to the working folder

        if [[ ! $line == make* ]]
        then
            cp $SWIS_path/bin/main.exe  $SWcode_executable
            echo Executable moved to $SWcode_executable
            compiled=true
        fi
    fi

    if [[ $rule =~ "clean" ]]; then 
        exit
    fi

    #--- Executing ------------------------------------------------------

    echo -e '\n******* Spin-wave code ***********'
    # basisname=$(grep "basisname" inputcard_$XX.inp | awk '{print $3}' | cut -d '"' -f 2 | cut -d ',' -f 2)
    # pairfile=$(grep "pairfile" inputcard_$XX.inp | awk '{print $3}' | cut -d '"' -f 2 | cut -d ',' -f 2)
    basisname=$(grep "basisname" inputcard_$XX.inp | awk '{print $3}' | cut -d '"' -f 2 | cut -d ',' -f 2 | cut -d '=' -f 2)
    pairfile=$(grep "pairfile" inputcard_$XX.inp | awk '{print $3}' | cut -d '"' -f 2 | cut -d ',' -f 2 | cut -d '=' -f 2)
    spirit_input=$(grep "spirit_input" inputcard_$XX.inp | awk '{print $3}' | cut -d '"' -f 2 | cut -d ',' -f 2 | cut -d '=' -f 2)
    latticefile=$(grep "latticefile" inputcard_$XX.inp | awk '{print $3}' | cut -d '"' -f 2 | cut -d ',' -f 2 | cut -d '=' -f 2)

    runsuccess=true
    if [ $executing -eq 1 ]; then
        should_execute=false
        if [ $forceexecuting -eq 1 ] ||  $compiled ; then
            should_execute=true
        fi

        declare -a files_to_check=(
                        $SWcode_executable 
                        "inputcard_$XX.inp"
                        $basisname
                        $pairfile
                        $spirit_input
                        $latticefile
                        )
        # now loop through the above array
        for file in "${files_to_check[@]}"
        do
           if ! cmp --silent $file temp/$file
           then
               echo "Modified file:" $file
               should_execute=true
           fi
        done

        # Executing the program
        if $should_execute
        then
            echo -e "\nSpinwave dispersion program starting:"
            ./$SWcode_executable  inputcard_$XX\.inp

            #if the program has run successfully, make the copying and backups 
            if [ "$?" -eq "0" ]; then
                mv precession_$XX.dat disp_analy_$XX.dat disp_unfol_$XX.dat  dispersion_$XX.dat latticExt_$XX.dat kpath_$XX.dat ine_intensities_$XX.dat outputfiles/

                #Copying files to be check for modifications on the next run
                for file in "${files_to_check[@]}"
                do
                    echo $file
                    #Copy the file and the folder structuring
                    rsync -R $file temp/
                done
            else
                runsuccess=false
            fi
        else
            echo "Calculated data are up to date!"
        fi

    fi # executing == 1

#--- Plottings ---------------------------------------------------------
    #--- Lattice Plotting -----------------------------------------------

    if [ $lattice -eq 1 ]; then
        # gr lattice_sky.py inputfiles/spinconfig_$XX.txt inputfiles/lattice_$XX.dat spinconfig_$XX.png
        
        latticefile=$(grep "latticefile" inputcard_$XX.inp | awk '{print $3}' | cut -d '"' -f 2)
        outputpng="$(basename "$basisname" .ovf).png"
        # outputpng="$(basename "$basisname" .txt).png"
        # outputpng="$(basename "$basisname" .txt).html"

        # gr lattice_sky.py $basisname $latticefile $outputpng
        python lattice_sky.py "outputfiles/latticExt_$XX.dat" "" $outputpng
        
        if [ "$?" -eq "0" ]; then
            if [ $host == 'theospc47' ] ; then
                # xdg-open $outputpng
                :
            else
                open $outputpng
            fi
        fi
        # open lattice_sky.html
    fi

    #--- Occupation number ----------------------------------------------

    if [ $occupation -eq 1 ]; then
        gnuplot occupation.gnu
        open occupation.html
    fi	

    #--- Dispersion -----------------------------------------------------
    if [ $dispersionplot -eq 1 ] && $runsuccess ; then
        echo
        echo "Plotting dispersion ..."

        # gnuscript="disp_unfol_multiPNG4panelsTHIN.gnu"


        maxomega=$(grep "maxomega" inputcard_$XX.inp | awk '{print $3}' | cut -d ',' -f 1)
        minomega=$(grep "minomega" inputcard_$XX.inp | awk '{print $6}' | cut -d ',' -f 1)

        gnuplot -e "specifier='$XX'; ymax='$maxomega'; ymin='$minomega'" $gnuscript

        #If the last command was sucessiful
        if [ "$?" -eq "0" ]; then

            if [ $host == 'theospc47' ] ; then
                xdg-open disp_unfolp1_$XX.png
                xdg-open disp_unfolp2_$XX.png
                xdg-open disp_unfolp3_$XX.png
                # :
            else
                open disp_unfolp1_$XX.png
                open disp_unfolp2_$XX.png
                open disp_unfolp3_$XX.png
            fi

        fi
            
        echo "Dispersion plotted."
    fi   

    #--- kpath -----------------------------------------------------
    if [ $kpath -eq 1 ] && $runsuccess ; then
        python kpath.py outputfiles/kpath_$XX.dat $XX
        
        if [ $host == 'theospc47' ] ; then
            xdg-open kpath_$XX.png
        else
            open kpath_$XX.png
        fi
        
        echo
        echo "Kpath plotted."
    fi
