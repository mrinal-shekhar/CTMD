#!/bin/bash
export AMBERHOME=/programs/x86_64-linux/system/sbgrid_bin
export input_folder="amber_MD_template"
#export CUDA_VISIBLE_DEVICES="0"

# SYSTEM
export prmtop=system.prmtop
export rd="." #/data/mshekhar/ETB/Global/MD_Amber/8XVH_gpcr_oxidized/MD
run_num=07
pre_num=06p5

NUM_REPEATS=3

rm -f COLVAR HILLS 07_Prod.* # No continuation
for i in `seq 1 $NUM_REPEATS`
do
	if test -f "COLVAR_$i"
	then
		echo "COLVAR exists for iteration $i. Skipping"
		continue
	fi

	pmemd.cuda \
-O \
-i $input_folder/07_Prod.in \
-o $rd/${run_num}_Prod.out \
-p $prmtop \
-c $rd/${pre_num}_Prod.rst \
-r $rd/${run_num}_Prod.rst \
-x $rd/${run_num}_Prod.nc \
-inf $rd/${run_num}_Prod.mdinfo

	mv COLVAR COLVAR_"$i"
	mv HILLS HILLS_"$i"
	mv 07_Prod.nc 07_Prod_iter"$i".nc
	mv 07_Prod.rst 07_Prod_iter"$i".rst
	mv 07_Prod.out 07_Prod_iter"$i".out
	rm -f 07_Prod.* 
done
