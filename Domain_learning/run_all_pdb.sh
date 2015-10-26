#!/bin/bash
current="HTH3"

make
if ! [ -d "all_vs_all/$current" ]
then
	mkdir all_vs_all/$current
fi
touch $current.txt
touch $current\_output.txt

date +%T

pdbs=(`ls ../restrict_domains/$current\_new`)
count=${#pdbs[@]}

for ((i=0; i<=count-1; i++ ))
do
	for ((j=i; j<= count-1; j++ ))
	do
		./align ../restrict_domains/$current\_new/${pdbs[i]} ../restrict_domains/$current\_new/${pdbs[j]} all_vs_all/$current/${pdbs[i]/.pdb/}_${pdbs[j]} $current
	done
done

date +%T
