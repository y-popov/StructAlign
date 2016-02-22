
from argparse import ArgumentParser
from subprocess import call
from sys import stdin
from os import system, access, F_OK, remove, devnull
from random import choice
from string import ascii_uppercase, digits



parser=ArgumentParser()
parser.add_argument('-i1', '--input1', help='First pdb-file')
parser.add_argument('-i2', '--input2', help='Second pdb-file')
parser.add_argument('-o', '--output', help='Name of output files, pdb1_pdb2 by default')
parser.add_argument('-c1', '--chain1', default='@', help='Protein chain in first pdb-file')
parser.add_argument('-c2', '--chain2', default='@', help='Protein chain in second pdb-file')
parser.add_argument('-s', '--supress', action='store_true', help='Supress program output text')

options=parser.parse_args()

chain1 = options.chain1.upper()
chain2 = options.chain2.upper()

if not options.output:
	output = options.input1[-8:-4]+'*_'+options.input2[-8:-4]+'*'
else:
        output = options.output

if not access(options.input1, F_OK):
	print ( "PDB entry <B>"+options.input1+"</B> does not exist" )
	exit(1)
if not access(options.input2, F_OK):
	print ( "PDB entry <B>"+options.input2+"</B> does not exist" )
	exit(1)

random_name = ''.join(choice(ascii_uppercase + digits) for i in range(10))
open(random_name+'.txt', 'w').close()

#devnull = open(devnull, 'w')
#args = 'algorithm.exe {} {} {}.pdb {} {} {}.txt'.format(options.input1, options.input2, output, chain1, chain2, random_name)
if options.supress:
        system('./align {} {} {}.pdb {} {} {}.txt > /dev/null'.format(options.input1, options.input2, output, chain1, chain2, random_name))
        #call(args, shell=False, stdout=devnull)
else:
        system('./align {} {} {}.pdb {} {} {}.txt'.format(options.input1, options.input2, output, chain1, chain2, random_name))
        #call(args, shell=False)

				
								
try:
	max_score = open("{}.txt".format(random_name), 'r')
except IOError as e:
	print "Sorry, the program has fault"
	print "I/O error({0}): {1}".format(e.errno, e.strerror)

				
max_score = max_score.read().splitlines()
for index, line in enumerate(max_score):
	if line.startswith("Error"):
		print "Error"
		print max_score[index+1]

	else:
		chain1, chain2, dna_chainA1, dna_chainA2, dna_chainB1, dna_chainB2, startA1, endA1, startA2, endA2, startB1, endB1, startB2, endB2, maxA, maxB, maxAc, maxBc = max_score[0].split()
		score = float( max_score[1] )
		dna1_chain1, dna1_chain2 = max_score[2], max_score[3] #sequence
		dna2_chain1, dna2_chain2 = max_score[4], max_score[5] #sequence
		dna11 = range(int(startA1), int(endA1)+1)
		dna12 = range(int(endA2), int(startA2)-1, -1)[::-1]
		dna21 = range(int(startB1), int(endB1)+1)
		dna22 = range(int(endB2), int(startB2)-1, -1)[::-1]

		code1 = options.input1[0:4]
		code2 = options.input2[0:4]

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
		

		print "Structures {} chain {} and {} chain {} were aligned with score: {}".format(code1, chain1, code2, chain2, score)
		print "The nucleotides with max measure: {}.{}:{}, {}.{}:{}".format(code1, dna_chainA1, maxA, code2, dna_chainB1, maxB)
		
		
		fasta = open(output[:-7]+output[-6:-1]+'.fasta', 'w')
		fasta.write(">{}.{} {}-{}\n{}\n\n>{}.{} {}-{}\n{}\n\n>{}.{} {}-{}\n{}\n\n>{}.{} {}-{}\n{}\n".format(code1, dna_chainA1, startA1, endA1, dna1_chain1, code2, dna_chainB1, startB1, endB1, dna2_chain1, code1, dna_chainA2, startA2, endA2, dna1_chain2, code2, dna_chainB2, startB2, endB2, dna2_chain2) )
		fasta.close()


		print 'To avoid possible overlapping of chain names, they were changed. The resulting pdb-file contains new names!'
		print "DNA chains:"
		print description_string

		print "Protein chains:"
		print "{}.{} -> E".format(code1, chain1)
		print "{}.{} -> F".format(code2, chain2)
		break

remove("{}.txt".format(random_name))
