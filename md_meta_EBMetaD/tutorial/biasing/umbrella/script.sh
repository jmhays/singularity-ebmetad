#!/bin/bash
GROMPP="/Users/chicco/Programs/gronamd/md_meta_devel/codes/gromacs4/gromacs-4.0.7/install_dir/bin/grompp_d"
MDRUN="/Users/chicco/Programs/gronamd/md_meta_devel/codes/gromacs4/gromacs-4.0.7/install_dir/bin/mdrun_d"
nstep=25
incr=0
ks=500
pi=3.14159265358979
start_config="2ala.gro"
input_file="md.mdp"
topol="gromacs.top"
indexfile="index.ndx"
###########################################################################
# here starts the umbrella sampling run 
###########################################################################
rm -rf metadatafile CV_*
for i in  `seq 0 $nstep`
    do
    val1=`echo $i\*2\* $pi/$nstep - $pi   | bc -l | awk '{printf("%12.6f",$1)}'`
    for j in  `seq 0 $nstep`
        do
        # 0 go forth 1 go back
        dir=`echo ${i}%2 | bc `
        if  [ "$dir" -eq "0" ]; then
        val2=`echo $j\*2\* $pi /$nstep - $pi | bc -l | awk '{printf("%12.6f",$1)}' `
        else 
        val2=`echo $pi -$j\*2\* $pi /$nstep  | bc -l | awk '{printf("%12.6f",$1)}'`
        fi
        incr=$((incr + 1)) 
        #
        # at the very beginning do first a steering where you intend to start the dynamics
        #
        if [ "$incr" -eq "1" ]; then
        cat >metadyn.dat <<EOF
        PRINT W_STRIDE 10 
        TORSION LIST 5 7 9 15   
        TORSION LIST 7 9 15 17 
        UMBRELLA CV 1 KAPPA ${ks} AT ${val1}
        UMBRELLA CV 2 KAPPA ${ks} AT ${val2}
        ENDMETA
EOF
        #
        # preproc 
        #
        $GROMPP -f md.mdp -c 2ala.gro -p gromacs.top -n index.ndx &>/dev/null 
        $MDRUN -plumed metadyn &>/dev/null
        #exit
        fi
        echo "ANG1 $val1 ANG2 $val2 INCR $incr"
        cat >metadyn.dat <<EOF
        PRINT W_STRIDE 10 
        TORSION LIST 5 7 9 15   
        TORSION LIST 7 9 15 17 
        UMBRELLA CV 1 KAPPA ${ks} AT ${val1}
        UMBRELLA CV 2 KAPPA ${ks} AT ${val2}
        ENDMETA
EOF
        $GROMPP -f md.mdp -c confout.gro -p gromacs.top -n index.ndx &>/dev/null 
        $MDRUN -plumed metadyn &>/dev/null
        grep -v FIELDS COLVAR | awk '{printf("%12.6f %12.6f %12.6f\n",$1,$2,$3)}'  >CV_${incr}
        echo "CV_${incr} ${val1} ${val2} ${ks} ${ks} ">>metadatafile 
        rm -rf \#*
        done
done
