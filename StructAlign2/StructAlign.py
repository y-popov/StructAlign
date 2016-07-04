
from argparse import ArgumentParser
from subprocess import check_output
#from sys import stdin
from os import system, access, F_OK, remove, devnull, path, makedirs, rename
from sys import argv
from random import choice
from string import ascii_uppercase, digits
from urllib2 import urlopen



parser=ArgumentParser(description="Align two DNA-protein complexes")
parser.add_argument('pdb1', help='First pdb-file')
parser.add_argument('pdb2', help='Second pdb-file')
parser.add_argument('-o', '--output', help='Output directory')
parser.add_argument('-c1', '--chain1', default='@', help='Protein chain in first pdb-file ')
parser.add_argument('-c2', '--chain2', default='@', help='Protein chain in second pdb-file')
parser.add_argument('-r1', '--range1', default='-', help='Use this range of protein in first pdb-file. Example: 5-60, -60, 5-')
parser.add_argument('-r2', '--range2', default='-', help='Use this range of protein in second pdb-file')
parser.add_argument('-i', '--internal', action="store_true", help="Use internal algorithm for complement nucleotide search; not stable but doesn't require 3DNA tools. DOES NOT WORK NOW.")
parser.add_argument('-s', '--supress', action='store_true', help='Suppress program internal output text')
parser.add_argument('-ss', '--supressAll', action='store_true', help='Suppress all program output text')

options=parser.parse_args()

not_installed = int( check_output(['{}./check_3dna.sh'.format(argv[0].replace("StructAlign.py", ''))]) )
if (not_installed==1) and (not options.internal):
	print "It seems you don't have installed 3DNA package! Please, install it or use '-i' option."
	print "You can download 3DNA from http://forum.x3dna.org/downloads/3dna-download/"
	print "Follow installation instructions from http://forum.x3dna.org/howtos/how-to-install-3dna-on-linux-and-windows/"
	exit(1)


chain1 = options.chain1.upper()
chain2 = options.chain2.upper()

def readRange(r):
	if '-' not in r:
		print "There is an error in range input. You missed '-'. Aborting..."
		exit(1)
	buff = r.split('-')
	if buff[0] == '':
		buff[0] = "zero"
	else:
		if not buff[0].lstrip('-').isdigit():
			print "The start of range is not numeric! Try again. Aborting.."
			exit(1)
	if buff[1] == '':
		buff[1] = "inf"
	else:
		if not buff[1].lstrip('-').isdigit():
			print "The end of range is not numeric! Try again. Aborting.."
			exit(1)
	return buff

range1 = readRange(options.range1)
range2 = readRange(options.range2)

code1 = options.pdb1[options.pdb1.rfind('/')+1:-4]
code2 = options.pdb2[options.pdb2.rfind('/')+1:-4]

if not options.output:
	output = code1+'@_'+code2+'@'
else:
	if not path.exists(options.output):
		makedirs(options.output)
        output = options.output.rstrip('/')+'/'+code1+'@_'+code2+'@'

dev_null = ''
dev_null_err = ''
if options.supressAll:
	options.supress = True
if options.supress:
	dev_null = ' > /dev/null'
	dev_null_err = ' 2&> /dev/null'

if not access(options.pdb1, F_OK):
	if not options.supressAll:
		print 'You do not have %s pdb-file! Downloading it...' % code1
	try:
		response = urlopen("http://files.rcsb.org/download/{}.pdb".format(code1))
		with open(options.pdb1, 'w') as dl:
			dl.write(response.read())
	except Exception:
		print "... aborting. PDB entry "+code1+" does not exist or you have not Internet connection!"
		exit(1)

if not access(options.pdb2, F_OK):
	if not options.supressAll:
		print 'You do not have %s pdb-file! Downloading it...' % code2
	try:
		response = urlopen("http://files.rcsb.org/download/{}.pdb".format(code2))
		with open(options.pdb2, 'w') as dl:
			dl.write(response.read())
	except Exception:
		print "... aborting. PDB entry "+code2+" does not exist or you have not Internet connection!"
		exit(1)

if path.isfile("out"):
	rename("out", "out.backup")
	if not options.supressAll:
		print "File 'out' was renamed to 'out.backup'!"
if path.exists("out"):
	print "Please, remove 'out' folder. It conflicts with program!"
	print 'Aborting'
	exit(1)

random_name = ''.join(choice(ascii_uppercase + digits) for i in range(10))
open(random_name+'.txt', 'w').close()

#devnull = open(devnull, 'w')
#args = 'algorithm.exe {} {} {}.pdb {} {} {}.txt'.format(options.pdb1, options.pdb2, output, chain1, chain2, random_name)

#print '{}./align {} {} {}.pdb {} {} {} {} {} {} {}.txt 0{}'.format(argv[0].replace("StructAlign.py", ''), options.pdb1, options.pdb2, output, chain1, chain2, range1[0], range1[1], range2[0], range2[1], random_name, dev_null)
system('{}./align {} {} {}.pdb {} {} {} {} {} {} {}.txt 0{}'.format(argv[0].replace("StructAlign.py", ''), options.pdb1, options.pdb2, output, chain1, chain2, range1[0], range1[1], range2[0], range2[1], random_name, dev_null))
        
#0 is for is_server
				
								
try:
	max_score = open("{}.txt".format(random_name), 'r')
	max_score = max_score.read().splitlines()
	if not max_score:
		print "Sorry, the program has fault"
		remove("{}.txt".format(random_name))
		exit(1)
except IOError as e:
	print "Sorry, the program has fault"
	print "I/O error({0}): {1}".format(e.errno, e.strerror)
	exit(1)


print ''
index = -1				
while index < len(max_score):
	index += 1
	line = max_score[index]
	if line.startswith("Error"):
		print "Error"
		print max_score[index+1]
		break
	elif line.startswith("Warning"):
		print line
		del max_score[index]
		index -= 1
	else:
		#print max_score[0]
		chain1, chain2, dna_chainA1, dna_chainA2, dna_chainB1, dna_chainB2, maxA, maxB, maxAc, maxBc, isreverse1, isreverse2 = max_score[0].split()
		score = float( max_score[1] )
		dna1_chain1, dna11 = max_score[2].split()[0], [int(x) for x  in max_score[2].split()[1].split(',')]
		dna1_chain2, dna12 = max_score[3].split()[0], [int(x) for x  in max_score[3].split()[1].split(',')]
		dna2_chain1, dna21 = max_score[4].split()[0], [int(x) for x  in max_score[4].split()[1].split(',')]
		dna2_chain2, dna22 = max_score[5].split()[0], [int(x) for x  in max_score[5].split()[1].split(',')] #sequence
		#dna11 = range(int(startA1), int(endA1)+1) if int(startA1)<int(endA1) else range(int(startA1), int(endA1)-1, -1)
		#dna12 = range(int(endA2), int(startA2)-1, -1)[::-1] if int(endA2)>int(startA2) else range(int(endA2), int(startA2)+1)[::-1]
		#dna21 = range(int(startB1), int(endB1)+1) if int(startB1)<int(endB1) else range(int(startB1), int(endB1)-1, -1)
		#dna22 = range(int(endB2), int(startB2)-1, -1)[::-1] if int(endB2)>int(startB2) else range(int(endB2), int(startB2)+1)[::-1]
		
		#print dna11, '\n', dna12, '\n', dna21, '\n', dna22, '\n'
		
		startA1, endA1 = dna11[0], dna11[-1]
		startA2, endA2 = dna12[0], dna12[-1]
		startB1, endB1 = dna21[0], dna21[-1]
		startB2, endB2 = dna22[0], dna22[-1]

		description_string = "{0}.{1} -> A\n{0}.{2} -> B\n{3}.{4} -> C\n{3}.{5} -> D".format(code1, dna_chainA1, dna_chainA2, code2, dna_chainB1, dna_chainB2)
		#if int(maxA) in dna12:
		if isreverse1 == '1':
			dna11, dna12 = dna12, dna11
			startA1, endA1, startA2, endA2 = startA2, endA2, startA1, endA1
			dna1_chain1, dna1_chain2 = dna1_chain2, dna1_chain1
			dna_chainA1, dna_chainA2 = dna_chainA2, dna_chainA1
			#maxA, maxAc = maxAc, maxA
		#if int(maxB) in dna22:
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
		
		if not options.supressAll:
			print "\nStructures {} chain {} and {} chain {} were aligned with score: {}".format(code1, chain1, code2, chain2, score)
			print "\nThe nucleotides with max measure: {}.{}:{}, {}.{}:{}".format(code1, dna_chainA1, maxA, code2, dna_chainB1, maxB)
		
		fasta_name = output.replace('@', '')+'.fasta'
		pdb_name = output[:output.find('@')]+chain1+output[output.find('@')+1:-1]+chain2+'.pdb'
		rename(output+'.pdb', pdb_name)
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
