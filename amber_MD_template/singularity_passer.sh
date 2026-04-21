#!/bin/bash

cd /scratch/$1
pwd
source /usr/software/sourceall.sh

if test $2 -eq 0
then
	# Setup MD
	./amber_MD_template/setup_MD.sh
else
	# Run BPMetad
	./amber_MD_template/run_prod.sh
fi
