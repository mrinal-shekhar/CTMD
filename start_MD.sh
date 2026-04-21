#!/bin/bash
#
fol=$1
if test -z "$fol"
then
	echo "Pick a folder name with files prepared"
	exit 1
fi

cd $fol
ln -s ../amber_MD_template/ .
./amber_MD_template/submit_MD.sh # Convert this to "sbatch ./amber_MD_template/submit_MD.sh" for clusters and ensure that the header is properly set

cd ..
echo "Completed"
