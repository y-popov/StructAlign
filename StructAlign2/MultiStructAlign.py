
from argparse import ArgumentParser
from subprocess import check_output
from sys import stdin
from os import system, access, F_OK, remove, devnull, path, makedirs, listdir, rename, getcwd, stat
from sys import argv
from random import choice
from string import ascii_uppercase, digits
from urllib2 import urlopen
from shutil import rmtree
from Bio import Phylo 
from Bio.Phylo import TreeConstruction
from pylab import show, savefig

def readRange(r, code, chain):
	#print code, chain, r
	if code not in ranges: 
		buff = ["zero", "inf"]
	elif chain not in ranges[code] and ranges[code].values() != ['@']:#len(ranges[code]) > 1:
		buff = ["zero", "inf"]
	else:
		if ranges[code].values() == ['@']:
			chain = ranges[code].keys()[0]
		buff = ranges[code][chain]
		if buff[0] == '':
			buff[0] = "zero"
		elif buff[0] != "zero":
			if not buff[0].lstrip('-').isdigit():
				print "The start of range is not numeric! Try again. Aborting.."
				exit(1)
		if buff[1] == '':
			buff[1] = "inf"
		elif buff[1] != "inf":
			if not buff[1].lstrip('-').isdigit():
				print "The end of range is not numeric! Try again. Aborting.."
				exit(1)
	#print code, chain, buff
	return buff

def StructAlign(pdb1, pdb2, outfile, warns, pairs, maxMs, ranges):
	code1 = pdb1[pdb1.rfind('/')+1:-1]
	code2 = pdb2[pdb2.rfind('/')+1:-1]
	
	chain1 = chain1_old = pdb1[-1].upper()
	chain2 = chain2_old = pdb2[-1].upper()
	
	range1 = readRange(ranges, code1, chain1)
	range2 = readRange(ranges, code2, chain2)
	
	pdb1_name = pdb1[:-1]+'.pdb'
	pdb2_name = pdb2[:-1]+'.pdb'
	
	output = code1+'@_'+code2+'@'
	
	try:
		system('{}./align {} {} All/{}.pdb {} {} {} {} {} {} {}.txt 0 > /dev/null'.format(argv[0].replace("MultiStructAlign.py", ''), pdb1_name, pdb2_name, output, chain1, chain2, range1[0], range1[1], range2[0], range2[1], outfile))
	except Exception as e:
		error.append(e)
	#print '{}./align {} {} All/{}.pdb {} {} {} {} {} {} {}.txt 0 > /dev/null'.format(argv[0].replace("MultiStructAlign.py", ''), pdb1_name, pdb2_name, output, chain1, chain2, range1[0], range1[1], range2[0], range2[1], outfile)
	
	max_score = open("{}.txt".format(outfile), 'r')
	max_score = max_score.read().splitlines()
	if not max_score:
		#print "Sorry, the program has fault"
		return None, None, None
		#exit(1)
	
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
			maxMs[pair] = [(code1, dna_chainA1, dna1_chain1, dna11, maxA, isreverse1), (code1, dna_chainA2, dna1_chain2, dna12, maxAc), (code2, dna_chainB1, dna2_chain1, dna21, maxB, isreverse2), (code2, dna_chainB2, dna2_chain2, dna22, maxBc)]
			#print maxMs[pair]
			
			if chain1 != chain1_old and code1 in ranges:
				ranges[code1][chain1] = ranges[code1][chain1_old]
			if chain2 != chain2_old and code2 in ranges:
				ranges[code2][chain2] = ranges[code2][chain2_old]
			
			break
				
	return score, chain1, chain2
	
def print_matrix(name, d, list):
	print ''
	print name+':'
	leng = len(list[0])

	for i in range(len(list)):
		print list[i][list[i].rfind('/')+1:],
		for j in range(len(list)):
			if i>=j:
				s = d[list[i]][list[j]]
				print "{:>7.0f} ".format(float(s)),
			else:
				print " "*8,
		print ""
	print ' '*7,
	for i in list:
		print "{:7} ".format(i[i.rfind('/')+1:]),
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
	#index = 1 if l[0]==rep else 0
	maxM1 = maxM[0:2] if index==1 else maxM[2:] #select maxM of repr
	maxM2 = maxM[0:2] if index==0 else maxM[2:]
	target = l[index]
	
	if maxM2[0][5] == '0':
		chain1, chain2 = maxM2[0][1], maxM2[1][1]
	else:
		chain1, chain2 = maxM2[1][1], maxM2[0][1]

	#print "{}./realign {}.pdb {} {} {} {} All/{}.pdb {} {} {} {} {}".format(argv[0].replace("MultiStructAlign.py", ''), repr_pdb[:-1], repr_pdb[-1], maxM1[0][1], maxM1[0][4], maxM1[1][1], pair, index, pdb_name, target[-1], chain1, chain2)
	system("{}./realign {}.pdb {} {} {} {} All/{}.pdb {} {} {} {} {}".format(argv[0].replace("MultiStructAlign.py", ''), repr_pdb[:-1], repr_pdb[-1], maxM1[0][1], maxM1[0][4], maxM1[1][1], pair, index, pdb_name, target[-1], chain1, chain2))
	
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
	index = 1 if pair[0]==rep else 0
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
		'''print '========================'
		for x in fastas:
			print x
		print "fasta: ", fasta.__repr__()
		print "refas: ", ref_fasta.__repr__()'''
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
parser.add_argument('-r', '--ranges', nargs='+', help='Ranges of proteins in format pdb:chain:start-end. Chain is not necessary (pdb::start--end).', default=[])
parser.add_argument('-o', '--output', default="multi", help="Name for output files")

options=parser.parse_args()

not_installed = int( check_output(['{}./check_3dna.sh'.format(argv[0].replace("MultiStructAlign.py", ''))]) )
if (not_installed==1) and (not options.internal):
	print "It seems you don't have installed 3DNA package! Please, install it."
	print "You can download 3DNA from http://forum.x3dna.org/downloads/3dna-download/"
	print "Follow installation instructions from http://forum.x3dna.org/howtos/how-to-install-3dna-on-linux-and-windows/"
	exit(1)

temp = options.output[:options.output.rfind('/')+1]
if temp != "":
	if not path.exists(temp):
		makedirs(temp)
temp = options.output[options.output.rfind('/')+1:]
if temp == "":
	options.output += "multi"
multy_pdb = options.output+".pdb"
multy_fasta1 = options.output+"1.fasta"
multy_fasta2 = options.output+"2.fasta"
multy_png = options.output+".png"
multy_descr = options.output+".txt"

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
		elif chain not in chains[pdb]:
			chains[pdb].append(chain)

ranges = {}
for handle in options.ranges:
	pdb, chain, r = handle.split(':')
	chain = chain.upper()
	if chain == '':
		chain = '@'
	if '--' not in r:
		print "There is an error in range input. You missed '--'. Aborting..."
		exit(1)
	r = r.split('--')
	if pdb not in ranges:
		ranges[pdb] = {chain:r}
	else:
		ranges[pdb][chain] = r
			
pdbs = []
used = []
files = []
if options.folder:
	if path.exists(options.folder):
		files = [options.folder.rstrip('/')+'/'+x for x in listdir(options.folder)]
	else:
		print "There is no folder %s! Aborting..." % options.folder
		exit(1)	

for pdb in options.pdbs+files:
	code = pdb[::-1].replace(".pdb"[::-1], '', 1)[::-1]
	code_id = code[code.rfind('/')+1:]
	if code_id not in used:
		buff = code
		if not access(pdb, F_OK):
			print 'You do not have %s pdb-file! Downloading it...' % code
			try:
				response = urlopen("http://files.rcsb.org/download/{}.pdb".format(code[code.rfind('/')+1:]))
				if options.output[:options.output.rfind('/')+1]:
					buff = options.output[:options.output.rfind('/')+1]+code[code.rfind('/')+1:]
				
				if not path.exists(buff[:buff.rfind('/')+1]):
					makedirs(buff[:buff.rfind('/')+1])
				with open("{}.pdb".format(buff), 'w') as dl:
					dl.write(response.read())
				if stat("{}.pdb".format(buff)).st_size == 0:
					print "Structure {} downloaded with error! Please, try again or download it manually.".format(code[code.rfind('/')+1:])
					exit(1)
			except Exception:
				print "... aborting. PDB entry "+code+" does not exist or you have not Internet connection!"
				exit(1)
		code_id = buff[buff.rfind('/')+1:]
		if code_id in chains:
			for chain in chains[code_id]:
				pdbs.append(buff+chain)
		else:
			pdbs.append(buff+'@')
		used.append(code_id)

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
errors = []
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
		#print pdbs[i], pdbs[j], '--{:->3.2%}-->'.format( s/total ), 
		failure_c = 0
		while failure_c < 3:
			score, chain1, chain2 = StructAlign(pdbs[i], pdbs[j], random_name, warnings, pairs, maxMs, ranges)
			if score is not None:
				failure_c = 3
			else:
				failure_c += 1
				print '{:->3.2%}\t{} {}\t{}'.format( s/total, pdbs[i], pdbs[j], "fail")
		if score is not None:
			pdbs[i] = pdbs[i][:-1]+chain1
			pdbs[j] = pdbs[j][:-1]+chain2
			print '{:->3.2%}\t{} {}\t{}'.format( s/total, pdbs[i], pdbs[j], "done")
			if pdbs[i] not in scores:
				scores[pdbs[i]] = {}
				PhyloM.append([])
				scoresM.append([])
			scores[pdbs[i]][pdbs[j]] = score
			scoresM[i].append(score)
		else:
			print '{:->3.2%}\t{} {}\t{}'.format( s/total, pdbs[i], pdbs[j], "error")
			if pdbs[i] not in scores:
				scores[pdbs[i]] = {}
				PhyloM.append([])
				scoresM.append([])
			scores[pdbs[i]][pdbs[j]] = 0
			scoresM[i].append(0)
		
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
def hide_inner(node):
	if node.name.startswith("Inner"):
		return None
	else:
		return node.name
try:
	Phylo.draw(tree, label_func=hide_inner, do_show=False)
	#show()
	savefig(multy_png)

except:
	print "Error: cannot create a figure of the tree! There is no DISPLAY variable\n"

		
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


f = open(multy_descr, 'w')
f.write("Model	PDB_ID	Protein_Chain	DNA_Chain_1	DNA_Chain_2\n")
for index, i in enumerate(description.split('\n')[1:-1]):
	buff = i.split()
	f.write("{}\t{}\t{}\t{}\t{}\n".format(buff[1][:-1], buff[2][:-1].upper(), buff[2][-1], fasta_align1[index].chain, fasta_align2[index].chain))
f.close()

if errors:
	print '\n'.join(errors)

if warnings:
	print '\n'+warnings
print "\nFiles {}, {}, {}, {} and {} were generated".format(multy_pdb, multy_fasta1, multy_fasta2, multy_png, multy_descr)
print description

		
remove("{}.txt".format(random_name))
x3dna_files = ["bestpairs.pdb", "col_chains.scr", "hel_regions.pdb", "ref_frames.dat", "bp_order.dat", "col_helices.scr", "out"]
for file in x3dna_files:
	remove(file)	
#rmtree('All')

