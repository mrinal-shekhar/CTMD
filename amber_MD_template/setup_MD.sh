#!/bin/bash
export AMBERHOME=/programs/x86_64-linux/system/sbgrid_bin
export CUDA_VISIBLE_DEVICES="0"
export input_folder="amber_MD_template"

# SYSTEM
export prmtop=system.prmtop
export name=system

## Minimise ##

if test -f "01_Min.rst"
then
	echo "01_Min.rst exists. Skipping 01_Min"
else
	pmemd \
-O \
-i $input_folder/01_Min.in \
-o 01_Min.out \
-p $prmtop \
-c ${name}.inpcrd \
-r 01_Min.rst 
fi

if test -f "02_Min2.rst"
then
	echo "02_Min2.rst exists. Skipping 02_Min"
else
	pmemd.cuda \
-O \
-i $input_folder/02_Min2.in \
-o 02_Min2.out \
-p $prmtop \
-c 01_Min.rst \
-r 02_Min2.rst 
fi

### Equilibration ##

if test -f "03_Heat.rst"
then
	echo "03_Heat.rst exists. Skipping 03_Heat"
else
	pmemd.cuda \
-O \
-i $input_folder/03_Heat.in \
-o 03_Heat.out \
-p $prmtop \
-c 02_Min2.rst \
-r 03_Heat.rst \
-x 03_Heat.nc \
-ref 02_Min2.rst
fi

if test -f "04_Heat2.rst"
then
	echo "04_Heat2.rst exists. Skipping 04_Heat"
else
	pmemd.cuda \
-O \
-i $input_folder/04_Heat2.in \
-o 04_Heat2.out \
-p $prmtop \
-c 03_Heat.rst \
-r 04_Heat2.rst \
-x 04_Heat2.nc \
-ref 03_Heat.rst
fi

## Peptide backbone restrained

if test -f "05_Back.rst"
then
	echo "05_Back.rst exists. Skipping 05_Back"
else
	pmemd.cuda \
-O \
-i $input_folder/05_Back.in \
-o 05_Back.out \
-p $prmtop \
-c 04_Heat2.rst \
-r 05_Back.rst \
-x 05_Back.nc \
-ref 04_Heat2.rst \
-inf 05_Back.mdinfo
fi

## Peptide C-alpha atoms only restrained

if test -f "06_Calpha.rst"
then
	echo "06_Calpha.rst exists. Skipping 06_Calpha"
else
	pmemd.cuda \
-O \
-i $input_folder/06_Calpha.in \
-o 06_Calpha.out \
-p $prmtop \
-c 05_Back.rst \
-r 06_Calpha.rst \
-x 06_Calpha.nc \
-ref 05_Back.rst \
-inf 06_Calpha.mdinfo
fi

echo "Setup complete"
echo "Running equilibriation"

if test -f "06p5_Prod.rst"
then
	echo "06p5_Prod.rst exists. Skipped equilibriation"
else
pmemd.cuda \
-O \
-i $input_folder/06p5_Prod.in \
-o 06p5_Prod.out \
-p $prmtop \
-c 06_Calpha.rst \
-r 06p5_Prod.rst \
-x 06p5_Prod.nc \
-inf 06p5_Prod.mdinfo
fi

