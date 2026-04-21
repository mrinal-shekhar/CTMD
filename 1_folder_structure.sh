#!/bin/bash
#

read -p "Copy files from 'extracted_poses/' (y/n)? " ync
yn=`echo $ync|cut -c1`
if [ "$yn" == "y" -o "$yn" == "Y" ]; then
	echo $yn "was YES"
	cp extracted_poses/*.mol2 .
	rm hit_data_mol2.spc -f
	cp extracted_poses/hit_data_mol2.spc .
else
	echo $yn "was NO"
fi

for f in `awk -e '{print $1}' hit_data_mol2.spc`
do
	echo $f
	folname=`echo $f|sed s/".mol2"//g`
	if test -d "$folname"
	then
		echo "Folder exists"
	else
		mkdir $folname
	fi
	cp $f $folname/
	cp $f $folname/ligand_raw.mol2
	grep $f -m1 hit_data_mol2.spc|awk -e '{print $3}' > $folname/Qnet
done
