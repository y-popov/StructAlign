
from argparse import ArgumentParser
from subprocess import check_output
from sys import stdin
from os import system, access, F_OK, remove, devnull, path, makedirs
from sys import argv
from random import choice
from string import ascii_uppercase, digits



parser=ArgumentParser(description="Align two DNA-protein complexes")
parser.add_argument('pdb1', help='First pdb-file')
parser.add_argument('pdb2', help='Second pdb-file')
parser.add_argument('-o', '--output', help='Name of output directory')
parser.add_argument('-c1', '--chain1', default='@', help='Protein chain in first pdb-file ')
parser.add_argument('-c2', '--chain2', default='@', help='Protein chain in second pdb-file')
parser.add_argument('-i', '--internal', action="store_true", help="Use internal algorithm for complement nucleotide search; not stable but doesn't require 3DNA tools")
parser.add_argument('-s', '--supress', action='store_true', help='Suppress program internal output text')
parser.add_argument('-ss', '--supressAll', action='store_true', help='Suppress all program output text')

options=parser.parse_args()

not_installed = int( check_output(['./check_3dna.sh']) )
if (not_installed==1) and (not options.internal):
	print "It seems you don't have installed 3DNA package! Please, install it or use '-i' option."
	print "You can download 3DNA from http://forum.x3dna.org/downloads/3dna-download/"
	print "Follow installation instructions from http://forum.x3dna.org/howtos/how-to-install-3dna-on-linux-and-windows/"
	exit(1)


chain1 = options.chain1.upper()
chain2 = options.chain2.upper()

code1 = options.pdb1[-8:-4]
code2 = options.pdb2[-8:-4]

if not options.output:
	output = options.pdb1[-8:-4]+'*_'+options.pdb2[-8:-4]+'*'
else:
	if not path.exists(options.output):
		makedirs(options.output)
        output = options.output.rstrip('/')+'/'+options.pdb1[-8:-4]+'*_'+options.pdb2[-8:-4]+'*'

if not access(options.pdb1, F_OK):
	print ( "PDB entry "+code1+" does not exist" )
	exit(1)
if not access(options.pdb2, F_OK):
	print ( "PDB entry "+code2+" does not exist" )
	exit(1)

if options.supressAll:
	options.supress = True

random_name = ''.join(choice(ascii_uppercase + digits) for i in range(10))
open(random_name+'.txt', 'w').close()

#devnull = open(devnull, 'w')
#args = 'algorithm.exe {} {} {}.pdb {} {} {}.txt'.format(options.pdb1, options.pdb2, output, chain1, chain2, random_name)
if options.supress:
        system('{}./align {} {} {}.pdb {} {} {}.txt > /dev/null'.format(argv[0].rstrip("StructAlign.py"), options.pdb1, options.pdb2, output, chain1, chain2, random_name))
        
else:
        system('{}./align {} {} {}.pdb {} {} {}.txt'.format(argv[0].rstrip("StructAlign.py"), options.pdb1, options.pdb2, output, chain1, chain2, random_name))
        

				
								
try:
	max_score = open("{}.txt".format(random_name), 'r')
	max_score = max_score.read().splitlines()
	if not max_score:
		print "Sorry, the program has fault"
except IOError as e:
	print "Sorry, the program has fault"
	print "I/O error({0}): {1}".format(e.errno, e.strerror)

				

for index, line in enumerate(max_score):
	if line.startswith("Error"):
		print "Error"
		print max_score[index+1]
	elif line.startswith("Warning"):
		print '\n', line
		del max_score[index]
	else:
		chain1, chain2, dna_chainA1, dna_chainA2, dna_chainB1, dna_chainB2, startA1, endA1, startA2, endA2, startB1, endB1, startB2, endB2, maxA, maxB, maxAc, maxBc = max_score[0].split()
		score = float( max_score[1] )
		dna1_chain1, dna1_chain2 = max_score[2], max_score[3] #sequence
		dna2_chain1, dna2_chain2 = max_score[4], max_score[5] #sequence
		dna11 = range(int(startA1), int(endA1)+1)
		dna12 = range(int(endA2), int(startA2)-1, -1)[::-1]
		dna21 = range(int(startB1), int(endB1)+1)
		dna22 = range(int(endB2), int(startB2)-1, -1)[::-1]
		print dna11, '\n', dna12, '\n', dna21, '\n', dna22, '\n'

		description_string = "{0}.{1} -> A\n{0}.{2} -> B\n{3}.{4} -> C\n{3}.{5} -> D".format(code1, dna_chainA1, dna_chainA2, code2, dna_chainB1, dna_chainB2)
		if int(maxA) in dna12:
			dna11, dna12 = dna12, dna11
			startA1, endA1, startA2, endA2 = startA2, endA2, startA1, endA1
			dna1_chain1, dna1_chain2 = dna1_chain2, dna1_chain1
			dna_chainA1, dna_chainA2 = dna_chainA2, dna_chainA1
			#maxA, maxAc = maxAc, maxA
		if int(maxB) in dna22:
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
		
		if not options.supressAll:
			print "\nStructures {} chain {} and {} chain {} were aligned with score: {}".format(code1, chain1, code2, chain2, score)
			print "\nThe nucleotides with max measure: {}.{}:{}, {}.{}:{}".format(code1, dna_chainA1, maxA, code2, dna_chainB1, maxB)
		
		fasta_name = output[:-7]+output[-6:-1]+'.fasta'
		pdb_name = output[:-7]+chain1+output[-6:-1]+chain2+'.pdb'
		fasta = open(fasta_name, 'w')
		fasta.write(">{}.{} {}-{}\n{}\n\n>{}.{} {}-{}\n{}\n\n>{}.{} {}-{}\n{}\n\n>{}.{} {}-{}\n{}\n".format(code1, dna_chainA1, startA1, endA1, dna1_chain1, code2, dna_chainB1, startB1, endB1, dna2_chain1, code1, dna_chainA2, startA2, endA2, dna1_chain2, code2, dna_chainB2, startB2, endB2, dna2_chain2) )
		fasta.close()

		if not options.supressAll:
			print '\nTo avoid possible overlapping of chain names, they were changed. The resulting pdb-file contains new names!'
			print "DNA chains:"
			print description_string

			print "Protein chains:"
			print "{}.{} -> E".format(code1, chain1)
			print "{}.{} -> F".format(code2, chain2)
			print "\n{} and {} were generated".format(pdb_name, fasta_name)
		break

remove("{}.txt".format(random_name))
if not options.internal:
	x3dna_files = ["bestpairs.pdb", "col_chains.scr", "hel_regions.pdb", "ref_frames.dat", "bp_order.dat", "col_helices.scr", "out"]
	for file in x3dna_files:
		remove(file)
