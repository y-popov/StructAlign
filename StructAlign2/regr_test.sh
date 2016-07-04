make
python StructAlign.py pdb/1puf.pdb pdb/3hdd.pdb -ss -o test
python StructAlign.py pdb/3hdd.pdb pdb/1puf.pdb -ss -o test
python StructAlign.py pdb/1puf.pdb pdb/1puf.pdb -c1 b -ss -o test
python StructAlign.py pdb/1puf.pdb pdb/1puf.pdb -c2 b -ss -o test
python StructAlign.py pdb/2pe5.149392.2pe5.pdb1.pdb pdb/2pe5.pdb -ss -o test 
python StructAlign.py pdb/1awc.pdb pdb/1awc.pdb -c2 b -ss -o test 
python StructAlign.py pdb/1dh3.pdb pdb/1t2k.pdb -o test -ss

