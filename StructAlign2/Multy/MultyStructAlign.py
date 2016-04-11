
from argparse import ArgumentParser
from subprocess import check_output
from sys import stdin
from os import system, access, F_OK, remove, devnull, path, makedirs, listdir
from sys import argv
from random import choice
from string import ascii_uppercase, digits
from urllib2 import urlopen

def StructAlign(pdb1, pdb2, outfile, warns):
	code1 = pdb1[pdb1.rfind('/')+1:-1]
	code2 = pdb2[pdb2.rfind('/')+1:-1]
	
	chain1 = pdb1[-1].upper()
	chain2 = pdb2[-1].upper()
	
	pdb1_name = pdb1[:-1]+'.pdb'
	pdb2_name = pdb2[:-1]+'.pdb'
	
	system('{}./align {} {} {}.pdb {} {} {}.txt 0 > /dev/null'.format(argv[0].replace("MultyStructAlign.py", pdb1_name, pdb2_name, output, chain1, chain2, outfile))
	
	max_score = open("{}.txt".format(random_name), 'r')
	max_score = max_score.read().splitlines()
	if not max_score:
		print "Sorry, the program has fault"
		remove("{}.txt".format(random_name))
		exit(1)
	
	index = -1
	while index < len(max_score):
		index += 1
		line = max_score[index]
		if line.startswith("Error"):
			print "Error"
			print max_score[index+1]
			break
		elif line.startswith("Warning"):
			warns += line
			del max_score[index]
			index -= 1
		else:
			chain1, chain2, dna_chainA1, dna_chainA2, dna_chainB1, dna_chainB2, startA1, endA1, startA2, endA2, startB1, endB1, startB2, endB2, maxA, maxB, maxAc, maxBc, isreverse1, isreverse2 = max_score[0].split()
			score = float( max_score[1] )
			dna1_chain1, dna1_chain2 = max_score[2], max_score[3] #sequence
			dna2_chain1, dna2_chain2 = max_score[4], max_score[5] #sequence
			dna11 = range(int(startA1), int(endA1)+1)
			dna12 = range(int(endA2), int(startA2)-1, -1)[::-1]
			dna21 = range(int(startB1), int(endB1)+1)
			dna22 = range(int(endB2), int(startB2)-1, -1)[::-1]
	
	return
	



parser=ArgumentParser(description="Align several DNA-protein complexes")
parser.add_argument('pdbs', nargs='*', help='pdb-files to align')
parser.add_argument('-d', '--folder', help="Folder with pdb-files")
parser.add_argument('-c', '--chains', nargs='+', help='Protein chains in format pdb:chain. You can assign several chains for pdb (-c pdb1:chain1 pdb1:chain2)')

options=parser.parse_args()

not_installed = int( check_output(['{}../check_3dna.sh'.format(argv[0].replace("MultyStructAlign.py", ''))]) )
if (not_installed==1) and (not options.internal):
	print "It seems you don't have installed 3DNA package! Please, install it or use '-i' option."
	print "You can download 3DNA from http://forum.x3dna.org/downloads/3dna-download/"
	print "Follow installation instructions from http://forum.x3dna.org/howtos/how-to-install-3dna-on-linux-and-windows/"
	exit(1)

chains = {}
answer = ''
for pair in options.chains:
	pdb, chain = pair.split(':')
	if pdb+'.pdb' not in options.pdbs:
		print "\nYou didn't specify file {}.pdb!".format(pdb)
		while answer.lower() not in ['n', 'y', 'no', 'yes']:
			answer = raw_input("Do you wish to continue without it? [y/n/a] ")
		if answer.lower() in ['n', 'no']:
			print 'Aborting...'
			exit(1)
		if answer.lower() in ['a', 'add']:
			if pdb not in chains:
				chains[pdb] = [chain]
			else:
				chains[pdb].append(chain)
	else:
		if pdb not in chains:
			chains[pdb] = [chain]
		else:
			chains[pdb].append(chain)
			
pdbs = []
files = []
if options.folder:
	if path.exists(options.folder):
		files = [options.folder.rstrip('/')+'/'+x for x in listdir(options.folder)]
	else:
		print "There is no folder %s! Aborting..." % options.folder
		exit(1)	
for pdb in options.pdbs+files:
	code = pdb.replace('.pdb', '')
	if code in chains:
		for chain in chains[code]:
			pdbs.append(code+chain)
	else:
		pdbs.append(code+'@')
	if not access(pdb, F_OK):
		print 'You do not have %s pdb-file! Downloading it...' % code
		try:
			response = urlopen("http://files.rcsb.org/download/{}.pdb".format(code))
			with open("{}.pdb".format(code), 'w') as dl:
				dl.write(response.read())
		except Exception:
			print "... aborting. PDB entry "+code+" does not exist or you have not Internet connection!"
			exit(1)

if path.isfile("out"):
	rename("out", "out.backup")
	if not options.supressAll:
		print "File 'out' was renamed to 'out.backup'!"
if path.exists("out"):
	print "Please, remove 'out' folder. It conflicts with program!"
	print 'Aborting'
	exit(1)

leng = len(pdbs)
scores = {}
random_name = ''.join(choice(ascii_uppercase + digits) for i in range(10))
open(random_name+'.txt', 'w').close()

if not path.exists("All"):
	makedirs('All')
warnings = ''
for i in range(leng):
	for j in range(i+1, leng):
		print pdbs[i], pdbs[j]
		score = StructAlign(pdbs[i], pdbs[j], random_name, warnings)
		
		
	
		

