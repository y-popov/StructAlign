StructAlign
======================================================================================
1) Place to desired folder
2) Unzip (tar -xf StructAlign.tar)
3) type 'make'
4) run as 'python StructAlign.py pdb1 pdb2'
5) для работы необходимы 3DNA и пакет питона BioPhylo

Type 'python StructAlign.py -h' for more options


pdb1*_pdb2*.pdb
	pdb-файл с совмещением структур, где * - имя белковой цепи

pdb1_pdb2.fasta
	выравнивание ДНК; совмещаются первые 2 и последние 2 последовательности
	
Чтобы использовать на kodomo, нужно указать путь к 3DNA:
export PATH=${PATH}:/home/preps/golovin/progs/X3DNA/bin
export X3DNA=/home/preps/golovin/progs/X3DNA



MultiStructAlign
======================================================================================
Запуск: python MultiStructAlign.py pdb1.pdb pdb2.pdb -d pdb_folder -c pdb2:chain1 pdb2:chain2 -r pdb1::10-50 pdb2:chain1:30-60 -o alignment
