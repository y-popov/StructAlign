StructAlign2

1) Unzip
2) Place to desired folder
3) type 'make'
4) run as 'python StructAlign.py pdb1 pdb2'

Type 'python StructAlign.py -h' for more options

Выходные файлы:
!Для формирования правильных имён имена pdb-файлов должны заканчиваться на pdb-код!

pdb1*_pdb2*.pdb
	pdb-файл с совмещением структур, где * - имя белковой цепи

pdb1_pdb2.fasta
	выравнивание ДНК; совмещаются первые 2 и последние 2 последовательности
