d=`date -I`
tar -cf StructAlign-$d.tar StructAlign.py MultyStructAlign.py algorithm.c realign.c pdb.c pdb.h nts.c changes.log check_3dna.sh Makefile README.txt 
