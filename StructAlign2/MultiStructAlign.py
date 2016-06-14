
from argparse import ArgumentParser
from subprocess import check_output
from sys import stdin
from os import system, access, F_OK, remove, devnull, path, makedirs, listdir, rename, getcwd
from sys import argv
from random import choice
from string import ascii_uppercase, digits
from urllib2 import urlopen
from shutil import rmtree
from Bio import Phylo 
from Bio.Phylo import TreeConstruction
from pylab import show, savefig

def StructAlign(pdb1, pdb2, outfile, warns, pairs, maxMs):
	code1 = pdb1[pdb1.rfind('/')+1:-1]
	code2 = pdb2[pdb2.rfind('/')+1:-1]
	
	chain1 = pdb1[-1].upper()
	chain2 = pdb2[-1].upper()
	
	pdb1_name = pdb1[:-1]+'.pdb'
	pdb2_name = pdb2[:-1]+'.pdb'
	
	output = code1+'@_'+code2+'@'
	
	system('{}./align {} {} All/{}.pdb {} {} {}.txt 0 > /dev/null'.format(argv[0].replace("MultiStructAlign.py", ''), pdb1_name, pdb2_name, output, chain1, chain2, outfile))
	
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
			chain1, chain2, dna_chainA1, dna_chainA2, dna_chainB1, dna_chainB2, maxA, maxB, maxAc, maxBc, isreverse1, isreverse2 = max_score[0].split()
			score = float( max_score[1] )
			dna1_chain1, dna11 = max_score[2].split()[0], [int(x) for x  in max_score[2].split()[1].split(',')]
			dna1_chain2, dna12 = max_score[3].split()[0], [int(x) for x  in max_score[3].split()[1].split(',')]
			dna2_chain1, dna21 = max_score[4].split()[0], [int(x) for x  in max_score[4].split()[1].split(',')]
			dna2_chain2, dna22 = max_score[5].split()[0], [int(x) for x  in max_score[5].split()[1].split(',')] #sequence
			#print chain1, chain2, dna_chainA1, dna_chainA2, dna_chainB1, dna_chainB2, maxA, maxB, maxAc, maxBc, isreverse1, isreverse2
			if isreverse1 == '1':
				#maxA, maxAc = maxAc, maxA
				dna11, dna12 = dna12, dna11
				dna1_chain1, dna1_chain2 = dna1_chain2, dna1_chain1
				dna_chainA1, dna_chainA2 = dna_chainA2, dna_chainA1
			if isreverse2 == '1':
				#maxB, maxBc = maxBc, maxB
				dna21, dna22 = dna22, dna21
				dna2_chain1, dna2_chain2 = dna2_chain2, dna2_chain1
				dna_chainB1, dna_chainB2 = dna_chainB2, dna_chainB1
			#print chain1, chain2, dna_chainA1, dna_chainA2, dna_chainB1, dna_chainB2, maxA, maxB, maxAc, maxBc, isreverse1, isreverse2
			pair = output[:output.find('@')]+chain1+output[output.find('@')+1:-1]+chain2

			pdb_name = 'All/'+pair+'.pdb'
			rename('All/'+output+'.pdb', pdb_name)
			
			pairs.append(pair)
			maxMs[pair] = [(code1, dna_chainA1, dna1_chain1, dna11, maxA), (code1, dna_chainA2, dna1_chain2, dna12, maxAc), (code2, dna_chainB1, dna2_chain1, dna21, maxB), (code2, dna_chainB2, dna2_chain2, dna22, maxBc)]
			break
				
	return score, chain1, chain2
	
def print_matrix(name, d, list):
	print ''
	print name+':'
	leng = len(list[0])
	print ' '*7,
	for i in list:
		print "{:7} ".format(i[i.rfind('/')+1:]),
	print ""

	for i in range(len(list)):
		print list[i][list[i].rfind('/')+1:],
		for j in range(len(list)):
			if i>=j:
				s = d[list[i]][list[j]]
				print "{:>7.0f} ".format(float(s)),
			else:
				print " "*8,
		print ""

def find_repr(m):
	print "\nSearching for the representative... ",
	mins = [min(x) for x in m]
	best = max(mins)
	i = mins.index(best)
	
	print m.names[i]
	return m.names[i], i

def reAlign(rep, repr_pdb, pair, maxM, pdb_name, c):
	l = pair.split('_')
	index = 0 if l[1]==rep else 1
	maxM1 = maxM[0:2] if index==1 else maxM[2:] #select maxM of repr
	maxM2 = maxM[0:2] if index==0 else maxM[2:]
	target = l[index]
	
	#print "{}./realign {}.pdb {} {} {} {} All/{}.pdb {} {} {} {} {} > /dev/null".format(argv[0].replace("MultiStructAlign.py", ''), repr_pdb[:-1], repr_pdb[-1], maxM1[0][1], maxM1[0][4], maxM1[1][1], pair, index, pdb_name, target[-1], maxM2[0][1], maxM2[1][1])
	system("{}./realign {}.pdb {} {} {} {} All/{}.pdb {} {} {} {} {} > /dev/null".format(argv[0].replace("MultiStructAlign.py", ''), repr_pdb[:-1], repr_pdb[-1], maxM1[0][1], maxM1[0][4], maxM1[1][1], pair, index, pdb_name, target[-1], maxM2[0][1], maxM2[1][1]))
	
	return "Model {}: {}\n".format(c, target)

def writePDB(fl, c, rep, pair, repr_pdb, maxM, pdb_name):
	fl = open(pdb_name, 'a')
	fl.write("MODEL{:>9}\n".format(c))
	fl.close()
	descr = reAlign(rep, repr_pdb, pair, maxM, pdb_name, c)
	fl = open(pdb_name, 'a')
	fl.write("ENDMDL\n")
	fl.close()
	return descr

class Fasta(object):
	def __init__(self, pdb, chain, seq, num_seq, max_m, max_ref):
		self.pdb = pdb
		self.chain = chain
		self.seq = seq
		self.num_seq = num_seq
		self.max_m = int(max_m)
		self.max_ref = int(max_ref)
		self.gaps = [0, 0]
	def __str__(self):
		return ">{}.{} {}-{}\n{}\n".format(self.pdb, self.chain, self.num_seq[0], self.num_seq[-1], '-'*self.gaps[0]+self.seq+'-'*self.gaps[1])
	def leng(self):
		return len(self.seq)+sum(self.gaps)
	def __repr__(self):
		return "pdb={} chain={} seq={} num_seq={} max_m={} max_ref={} gaps={}".format(self.pdb, self.chain, self.seq, self.num_seq, self.max_m, self.max_ref, self.gaps)
	
def createFASTA(fastas1, fastas2, maxMs_all, rep, pair):
	maxMs = maxMs_all[pair]
	maxM_ref = maxMs_all[rep+'_'+rep][:2]
	pair = pair.split('_')
	#print pair, rep
	index = 0 if pair[1]==rep else 1
	#print index
	target = pair[index]
	
	if index == 0:
		maxM = maxMs[:2]
		maxM2 = maxMs[2:]
	else:
		maxM = maxMs[2:]
		maxM2 = maxMs[:2]
	for ind, i in enumerate(maxM):
		#print i, maxM2[ind][4]
		buff = Fasta(i[0], i[1], i[2], i[3], i[4], maxM2[ind][4])
		if maxM2[ind][1] == maxM_ref[ind][1]:
			if ind == 0:
				fastas1.append(buff)
			if ind == 1:
				fastas2.append(buff)
		else:
			if ind == 1:
				fastas1.append(buff)
			if ind == 0:
				fastas2.append(buff)

def alignFASTA(fastas, ref_fasta):
	for index, fasta in enumerate(fastas):
		"""print '========================'
		for x in fastas:
			print x
		print "fasta: ", fasta.__repr__()
		print "refas: ", ref_fasta.__repr__()"""
		delta = fasta.num_seq.index(fasta.max_m)+fasta.gaps[0] - ref_fasta.num_seq.index(fasta.max_ref)-ref_fasta.gaps[0]
		if delta < 0:
			fasta.gaps[0] -= delta
		if delta > 0:
			ref_fasta.gaps[0] += delta
			for i in range(index):
				fastas[i].gaps[0] += delta
		delta = fasta.leng() - ref_fasta.leng()
		if delta < 0:
			fasta.gaps[1] -= delta
		else:
			ref_fasta.gaps[1] += delta
			for i in range(index):
				fastas[i].gaps[1] += delta		
	
def writeFASTA(f, fasta):
	for handle in fasta:
		f.write(handle.__str__())


parser=ArgumentParser(description="Align several DNA-protein complexes")
parser.add_argument('pdbs', nargs='*', help='pdb-files to align')
parser.add_argument('-d', '--folder', help="Folder with pdb-files")
parser.add_argument('-c', '--chains', nargs='+', help='Protein chains in format pdb:chain. You can assign several chains for pdb (-c pdb1:chain1 pdb1:chain2)', default=[])
parser.add_argument('-o', '--output', default="multi", help="Name for output files")

options=parser.parse_args()

not_installed = int( check_output(['{}./check_3dna.sh'.format(argv[0].replace("MultiStructAlign.py", ''))]) )
if (not_installed==1) and (not options.internal):
	print "It seems you don't have installed 3DNA package! Please, install it."
	print "You can download 3DNA from http://forum.x3dna.org/downloads/3dna-download/"
	print "Follow installation instructions from http://forum.x3dna.org/howtos/how-to-install-3dna-on-linux-and-windows/"
	exit(1)

multy_pdb = options.output+".pdb"
multy_fasta1 = options.output+"1.fasta"
multy_fasta2 = options.output+"2.fasta"
multy_png = options.output+".png"

chains = {}
for pair in options.chains:
	pdb, chain = pair.split(':')
	if pdb+'.pdb' not in [x[x.rfind('/')+1:] for x in options.pdbs]:
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
	code = pdb[::-1].replace(".pdb"[::-1], '', 1)[::-1]
	if code[code.rfind('/')+1:] in chains:
		for chain in chains[code[code.rfind('/')+1:]]:
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
maxMs = {}
PhyloM = []
scoresM = []
s = 0
total = (leng**2+leng)/2.0
print ''
for i in range(leng):
	for j in range(i+1): #or i+1?
		s += 1
		print pdbs[i], pdbs[j], '--{:->3.2%}-->'.format( s/total ), 
		score, chain1, chain2 = StructAlign(pdbs[i], pdbs[j], random_name, warnings, pairs, maxMs)
		pdbs[i] = pdbs[i][:-1]+chain1
		pdbs[j] = pdbs[j][:-1]+chain2
		print pdbs[i], pdbs[j]
		if pdbs[i] not in scores:
			scores[pdbs[i]] = {}
			PhyloM.append([])
			scoresM.append([])
		scores[pdbs[i]][pdbs[j]] = score
		scoresM[i].append(score)
		
print_matrix("Scores", scores, pdbs)
scoresM = TreeConstruction._Matrix([x[x.rfind('/')+1:] for x in pdbs], scoresM)

distances = {}
for i in range(leng):
	distances[pdbs[i]] = {}
	for j in range(i+1):
		distances[pdbs[i]][pdbs[j]] = (scores[pdbs[i]][pdbs[i]]+scores[pdbs[j]][pdbs[j]])/2.0 - scores[pdbs[i]][pdbs[j]]
		PhyloM[i].append(distances[pdbs[i]][pdbs[j]])
PhyloM = TreeConstruction._DistanceMatrix([x[x.rfind('/')+1:] for x in pdbs], PhyloM)
print_matrix("Distances", distances, pdbs)

tree = TreeConstruction.DistanceTreeConstructor().upgma(PhyloM)
Phylo.draw_ascii(tree)
tree.ladderize()
#Phylo.draw_graphviz(tree, node_size=0)
Phylo.draw(tree, show_confidence=False, do_show=False)
#show()
savefig(multy_png)

		
repres, repr_index = find_repr(scoresM)


f = open(multy_pdb, 'w')
f.write( "HEADER{:>60}\n".format("MultiStructAlign") )
f.write( "TITLE{:>60}\n".format("title") )
f.close()

c = 1
description = "\n"
for pair in pairs:
	if repres in pair:
		#reAlign(repres, pdbs[repr_index], pair, maxMs[pair], multy_pdb)
		description += writePDB(f, c, repres, pair, pdbs[repr_index], maxMs[pair], multy_pdb)
		c += 1

f = open(multy_pdb, 'a')
f.write( "END\n" )
f.close()


fasta_align1 = [] #array of fasta-classes
fasta_align2 = []

temp = maxMs[repres+'_'+repres][0]
ref_fasta1 = Fasta(temp[0], temp[1], temp[2], temp[3], temp[4], temp[4])
temp = maxMs[repres+'_'+repres][1]
ref_fasta2 = Fasta(temp[0], temp[1], temp[2], temp[3], temp[4], temp[4])

for pair in pairs:
	if repres in pair:
		createFASTA(fasta_align1, fasta_align2, maxMs, repres, pair)

g = open(multy_fasta1 ,'w')
alignFASTA(fasta_align1, ref_fasta1)
writeFASTA(g, fasta_align1)
g.close()

g = open(multy_fasta2 ,'w')
alignFASTA(fasta_align2, ref_fasta2)
writeFASTA(g, fasta_align2)
g.close()


if warnings:
	print '\n'+warnings
print "\nFiles {}, {}, {} and {} were generated".format(multy_pdb, multy_fasta1, multy_fasta2, multy_png)
print description
		
remove("{}.txt".format(random_name))
x3dna_files = ["bestpairs.pdb", "col_chains.scr", "hel_regions.pdb", "ref_frames.dat", "bp_order.dat", "col_helices.scr", "out"]
for file in x3dna_files:
	remove(file)	
#rmtree('All')

