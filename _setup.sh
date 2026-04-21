#!/bin/bash
#

ligname=$1
q_net=$2
add_ions=$3
if test -z "$ligname"
then
	echo "Please provide ligand mol2 file"
	exit 1
fi
if test -f "ligand.frcmod"
then
	echo "Ligand parameters exist. Skipping antechamber"
else
	if test -z "$q_net"
	then
		q_net=`grep $ligname ../hit_data_mol2.spc|awk -e '{print $3}'`
		if test $? -ne 0
		then
			echo "$ligname not found in charge database ../hit_data_mol2.spc"
			exit 1
		fi
	fi
	antechamber -i $ligname -fi mol2 -fo mol2 -o ligand.mol2 -nc $q_net -at gaff2 -c bcc -rn LIG
	if test $? -ne 0
	then
		echo "ERR: antechamber failed!"
		exit 2
	fi
	parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod
	rm -f sqm.* ANTECHAMBER*
fi

rm -f leap.log
tleap -f ../setup.leap
vol=`grep Volume leap.log |grep "A^3"|cut -d" " -f4`
echo "Final volume: $vol"
echo "Atom count raw: mul $vol 9.033e-05"
nat=`mul $vol 9.033e-05|cut -d"." -f1` #9.033e-05 atoms/A^3 is target conc
cp ../setup.leap ./setup.leap
if test -n "$add_ions"
then
	nat_cl=$nat
	nat_na=$nat
else
	nat_cl=0
	nat_na=0
fi

if test $q_net -ne 0
then
	echo "Ligand was not neutral. Adding ions"
	add_ions=True
	if test $q_net -gt 0
	then
		nat_cl=`expr $nat_cl "+" $q_net`
	else
		nat_na=`expr $nat_na "-" $q_net` # q_net is -ve
	fi
fi

if test -n "$add_ions"
then
	echo "Adding ions manually"
	echo $nat sodium and chloride ions to be added

	sed /"^ *quit"/d setup.leap > modsetup.leap
	echo "addionsrand solvcomplex Cl- $nat_cl" >> modsetup.leap
	echo "addionsrand solvcomplex Na+ $nat_na" >> modsetup.leap
	echo "savepdb solvcomplex system.pdb" >> modsetup.leap
	echo "saveamberparm solvcomplex system.prmtop system.inpcrd" >> modsetup.leap
	echo "quit" >> modsetup.leap

	rm -f system.prmtop system.rst7 system.pdb
	tleap -f modsetup.leap
fi
