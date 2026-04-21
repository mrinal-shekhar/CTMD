#!/bin/bash
export AMBERHOME=/programs/x86_64-linux/system/sbgrid_bin
export CUDA_VISIBLE_DEVICES="0"

# SYSTEM
export prmtop=System_Prep/Combine/system.prmtop
export name=System_Prep/Combine/system
export hmass=system_mol1.hmass.prmtop

## Minimise ##
pmemd \
-O \
-i 01_Min.in \
-o 01_Min.out \
-p ../$prmtop \
-c ../${name}.inpcrd \
-r 01_Min.rst 

