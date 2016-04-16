
from argparse import ArgumentParser
from subprocess import check_output
from sys import stdin
from os import system, access, F_OK, remove, devnull, path, makedirs, listdir, rename, getcwd
from sys import argv
from random import choice
from string import ascii_uppercase, digits
from urllib2 import urlopen
from shutil import rmtree

def StructAlign(pdb1, pdb2, outfile, warns, pairs):
	code1 = pdb1[pdb1.rfind('/')+1:-1]
	code2 = pdb2[pdb2.rfind('/')+1:-1]
	
	chain1 = pdb1[-1].upper()
	chain2 = pdb2[-1].upper()
	
	pdb1_name = pdb1[:-1]+'.pdb'
	pdb2_name = pdb2[:-1]+'.pdb'
	
	output = code1+'@_'+code2+'@'
	
	system('{}./align {} {} All/{}.pdb {} {} {}.txt 0 > /dev/null'.format(argv[0].replace("MultyStructAlign.py", ''), pdb1_name, pdb2_name, output, chain1, chain2, outfile))
	
	max_score = open("{}.txt".format(outfile), 'r')
	max_score = max_score.read().splitlines()
	if not max_score:
		print "Sorry, the program has fault"
		remove("{}.txt".format(outfile))
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
			warns += line+'\n'
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
			
			if isreverse1 == '1':
				dna11, dna12 = dna12, dna11
				startA1, endA1, startA2, endA2 = startA2, endA2, startA1, endA1
				dna1_chain1, dna1_chain2 = dna1_chain2, dna1_chain1
				dna_chainA1, dna_chainA2 = dna_chainA2, dna_chainA1
			if isreverse2 == '1':
				dna21, dna22 = dna22, dna21
				startB1, endB1, startB2, endB2 = startB2, endB2, startB1, endB1
				dna2_chain1, dna2_chain2 = dna2_chain2, dna2_chain1
				dna_chainB1, dna_chainB2 = dna_chainB2, dna_chainB1
			#create alignment 1
			delta = dna11.index(int(maxA)) - dna21.index(int(maxB))
			if delta < 0:
				dna1_chain1 = '-'*(-delta)+dna1_chain1
			if delta > 0:
				dna2_chain1 = '-'*delta+dna2_chain1
			delta = len(dna1_chain1)-len(dna2_chain1)
			if delta < 0:
				dna1_chain1 += '-'*(-delta)
			else:
				dna2_chain1 += '-'*delta
			#create alignment 2
			delta = dna12.index(int(maxAc)) - dna22.index(int(maxBc))
			if delta < 0:
				dna1_chain2 = '-'*(-delta)+dna1_chain2
			if delta > 0:
				dna2_chain2 = '-'*delta+dna2_chain2
			delta = len(dna1_chain2)-len(dna2_chain2)
			if delta < 0:
				dna1_chain2 += '-'*(-delta)
			else:
				dna2_chain2 += '-'*delta
			
			pair = output[:output.find('@')]+chain1+output[output.find('@')+1:-1]+chain2
			fasta_name = 'All/'+pair+'.fasta'
			pdb_name = 'All/'+pair+'.pdb'
			rename('All/'+output+'.pdb', pdb_name)
			fasta = open(fasta_name, 'w')
			fasta.write(">{}.{} {}-{}\n{}\n\n>{}.{} {}-{}\n{}\n\n>{}.{} {}-{}\n{}\n\n>{}.{} {}-{}\n{}\n".format(code1, dna_chainA1, startA1, endA1, dna1_chain1, code2, dna_chainB1, startB1, endB1, dna2_chain1, code1, dna_chainA2, startA2, endA2, dna1_chain2, code2, dna_chainB2, startB2, endB2, dna2_chain2) )
			fasta.close()
			
			pairs.append(pair)
			break
				
	return score, chain1, chain2
	
def print_matrix(d, list):
	print ''
	leng = len(list[0])
	print ' '*7,
	for i in list:
		print "{:7} ".format(i[i.rfind('/')+1:]),
	print ""

	for i in range(len(list)):
		print list[i][list[i].rfind('/')+1:],
		for j in range(len(list)):
			if j>=i:
				s = d[list[i]][list[j]]
				print "{:>7.0f} ".format(float(s)),
			else:
				print " "*8,
		print ""

def find_repr(d, pdbs):
	return pdbs[0][pdbs[0].rfind('/')+1:]

def writePDB(fl, c, rep, pdb):
	fl.write("MODEL{:>9}\n".format(c))
	
	pair = pdb.split('_')
	index = 0 if pair[1]==rep else 1
	target = pair[index]

	if index == 0:
		chains = ['A', 'B', 'E']
	else:
		chains = ['C', 'D', 'F']
	
	pdb_file = open("All/"+pdb+".pdb", 'r')
	for line in pdb_file:
		if len(line)>21:
			if line[21] in chains:
				fl.write(line)
	pdb_file.close()
	
	fl.write("ENDMDL\n")

def writeFASTA(fl, rep, pdb):
	pair = pdb.split('_')
	index = 0 if pair[1]==rep else 1
	target = pair[index]

	fasta_file = open("All/"+pdb+".fasta", 'r')
	flag = 0
	count = []
	for line in fasta_file:
		if flag == 1:
			fl.write(line+'\n')
			flag = 0
		if line[1:7] not in count:
			if target[:-1]+'.' in line:
				flag = 1
				fl.write(line)
				count.append(line[1:7])
	
	fasta_file.close()

parser=ArgumentParser(description="Align several DNA-protein complexes")
parser.add_argument('pdbs', nargs='*', help='pdb-files to align')
parser.add_argument('-d', '--folder', help="Folder with pdb-files")
parser.add_argument('-c', '--chains', nargs='+', help='Protein chains in format pdb:chain. You can assign several chains for pdb (-c pdb1:chain1 pdb1:chain2)', default=[])

options=parser.parse_args()

not_installed = int( check_output(['{}../check_3dna.sh'.format(argv[0].replace("MultyStructAlign.py", ''))]) )
if (not_installed==1) and (not options.internal):
	print "It seems you don't have installed 3DNA package! Please, install it."
	print "You can download 3DNA from http://forum.x3dna.org/downloads/3dna-download/"
	print "Follow installation instructions from http://forum.x3dna.org/howtos/how-to-install-3dna-on-linux-and-windows/"
	exit(1)

chains = {}
for pair in options.chains:
	pdb, chain = pair.split(':')
	if pdb+'.pdb' not in options.pdbs:
		answer = ''
		print "\nYou didn't specify file {}.pdb!".format(pdb)
		while answer.lower() not in ['n', 'y', 'no', 'yes', 'a', 'add']:
			answer = raw_input("Do you wish to continue without it? [y/n/a] ")
		if answer.lower() in ['n', 'no']:
			print 'Aborting...'
			exit(1)
		if answer.lower() in ['a', 'add']:
			options.pdbs.append(pdb+'.pdb')
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
pairs = []
s = 0
total = (leng**2+leng)/2.0
print ''
for i in range(leng):
	for j in range(i, leng): #or i+1?
		s += 1
		print pdbs[i], pdbs[j], '--{:->3.2%}-->'.format( s/total ), 
		score, chain1, chain2 = StructAlign(pdbs[i], pdbs[j], random_name, warnings, pairs)
		pdbs[i] = pdbs[i][:-1]+chain1
		pdbs[j] = pdbs[j][:-1]+chain2
		print pdbs[i], pdbs[j]
		if j == i:
			scores[pdbs[i]] = {}
		scores[pdbs[i]][pdbs[j]] = score
		
print_matrix(scores, pdbs)	


repres = find_repr(scores, pdbs)

f = open('multy.pdb', 'w')
f.write( "HEADER{:>60}\n".format("MultyStructAlign") )
f.write( "TITLE{:>60}\n".format("title") )

c = 1
for pair in pairs:
	if repres in pair:
		writePDB(f, c, repres, pair)
		c += 1

f.write( "END\n" )
f.close()


g = open('multy.fasta' ,'w')
for pair in pairs:
	if repres in pair:
		writeFASTA(g, repres, pair)
g.close()

if warnings:
	print '\n'+warnings
print "\nFiles multy.pdb and multy.fasta were generated"
		
remove("{}.txt".format(random_name))
x3dna_files = ["bestpairs.pdb", "col_chains.scr", "hel_regions.pdb", "ref_frames.dat", "bp_order.dat", "col_helices.scr", "out"]
for file in x3dna_files:
	remove(file)	
#rmtree('All')

