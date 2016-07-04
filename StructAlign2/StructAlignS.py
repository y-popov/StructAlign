#!/usr/bin/python

from sys import stdin
from os import system, access, F_OK, rename
from random import choice
from string import ascii_uppercase, digits


def complement(x):
	alph = 'ATGC-?'
	compl_alph = 'TACG--'
	s = ''
	for i in x:
		s += compl_alph[alph.index(i)]
	return s


print "Content-type: text/html\n\n"
print "<HTML>\n<HEAD>"
print '  <script type="text/javascript" src="/tools/jsmol/Jmol.js"></script> \n \
  <script type="text/javascript"> \n \
    jmolInitialize("/tools/jsmol/"); \n \
  </script> \n' 
print '''<TITLE>Result</TITLE>
	<link href='https://fonts.googleapis.com/css?family=Roboto:400,300,500,100' rel='stylesheet' type='text/css'>
	<link href="../StructAlign/StructAlign.css" rel="stylesheet" type="text/css">
</HEAD>\n'''
print "<BODY>"
print '\
  <script type="text/javascript">    \n \
    function jmMessage ( apName, Value )  \n \
    {   \n \
      var el = document.getElementById ( "jmConsole" ); \n \
      el.innerHTML = el.innerHTML + "<br>" + Value;   \n \
    }   \n \
  </script> \n '


print '''<div class="header">
			
	</div>

	<div class="nav">
		<div class="header">
			<h3>StructAlign</h3> 		
			<p>A program for alignment of structures of DNA-protein complexes</p>
		</div>

		 
		<ul>
			<a href="http://mouse.belozersky.msu.ru/tools/StructAlign.html"><li>Home</li></a>
			<a href="http://mouse.belozersky.msu.ru/tools/StructAlign/StructAlign_feedback.html"><li>Feedback</li></a>
			<a href="http://mouse.belozersky.msu.ru/tools/StructAlign/StructAlign_help.html"><li>Help</li></a>
			<a href="http://mouse.belozersky.msu.ru/tools/StructAlign/StructAlign_about.html"><li>About</li></a>
		</ul>
	</div>'''

a = stdin.readlines()

for b in a:
	c = b.split("&")
	c.sort()
	for d in c:
		(name,value) = d.split("=")
		if name == "chain1":
			if value == "":
				chain1 = "@"
			else:
				chain1 = value.upper()
		elif name == "chain2":
			if value == "":
				chain2 = "@"
			else:
				chain2 = value.upper()
		
		elif name == "start1":
			start1 = value if value != "" else "zero"
		elif name == "end1":
			end1 = value if value != "" else "inf"
		elif name == "start2":
			start2 = value if value != "" else "zero"
		elif name == "end2":
			end2 = value if value != "" else "inf"

		elif name == "zpdbcode1":
			code1 = value[0:4].lower()
			filename = "/mnt/databanks/pdb/"+value[1:3].lower()+"/pdb"+code1+".ent.gz"
			if access(filename, F_OK): # this PDB entry exists
				if not access("../tmp/"+code1+".pdb", F_OK):
					system("gunzip -c " + filename + "> ../tmp/"+code1+".pdb")
			else:
				print ( "<H3>Error</H3>" )
				print ( "<p>PDB entry <B>"+code1+"</B> does not exist</p>" )
		elif name == "zpdbcode2":
			code2 = value[0:4].lower()
			filename = "/mnt/databanks/pdb/"+value[1:3].lower()+"/pdb"+code2+".ent.gz"
			if access(filename, F_OK):
				if not access("../tmp/"+code2+".pdb", F_OK):
					system("gunzip -c " + filename + " > ../tmp/"+code2+".pdb")


				random_name = ''.join(choice(ascii_uppercase + digits) for i in range(10))
				system("mkdir ../tmp/StructAlign/{}".format(random_name))
				system("chmod a=rwx ../tmp/StructAlign/{}".format(random_name))
				system("touch ../tmp/StructAlign/{}/result.txt".format(random_name) )
				#system("StructAlign/3dna.sh 1puf.pdb %s" % random_name)
				#system("chmod a=rwx ../tmp/StructAlign/{}/result.txt".format(random_name) )
				
				#system("/home/popov/bin/align ../tmp/{0}.pdb ../tmp/{1}.pdb ../tmp/StructAlign/{0}*_{1}*.pdb {2} {3} {4} > ../tmp/StructAlign/log.txt".format(code1, code2, chain1, chain2, random_name) )
				system("StructAlign/align ../tmp/{0}.pdb ../tmp/{1}.pdb ../tmp/StructAlign/{8}/{0}@_{1}@.pdb {2} {3} {4} {5} {6} {7} /var/www/tools/tmp/StructAlign/{8}/result.txt 1 >> ../tmp/StructAlign/{8}/log.txt".format(code1, code2, chain1, chain2, start1, end1, start2, end2, random_name) ) #1 is for is_server				
				try:
					max_score = open("../tmp/StructAlign/{}/result.txt".format(random_name), 'r')
					max_score = max_score.read().splitlines()
					if not max_score:
						print "<p>Sorry, the program has fault</p>"
				except IOError as e:
					print "<p>Sorry, the program has fault</p>"
					print "I/O error({0}): {1}".format(e.errno, e.strerror)

				warnings = ''
				index = -1
				while index < len(max_score):
					index += 1
					line  = max_score[index]
					if line.startswith("Error"):
						print "<H3>Error</H3>"
						print "<p> "+max_score[index+1]+" </p>"
					elif line.startswith("Warning"):
						warnings += '<p>'+line+'</p>\n'
						del max_score[index]
						index -= 1

					else:
						chain1, chain2, dna_chainA1, dna_chainA2, dna_chainB1, dna_chainB2, maxA, maxB, maxAc, maxBc, isreverse1, isreverse2 = max_score[0].split()
						score = float( max_score[1] )
						dna1_chain1, dna11 = max_score[2].split()[0], [int(x) for x  in max_score[2].split()[1].split(',')]
						dna1_chain2, dna12 = max_score[3].split()[0], [int(x) for x  in max_score[3].split()[1].split(',')]
						dna2_chain1, dna21 = max_score[4].split()[0], [int(x) for x  in max_score[4].split()[1].split(',')]
						dna2_chain2, dna22 = max_score[5].split()[0], [int(x) for x  in max_score[5].split()[1].split(',')]

						startA1, endA1 = dna11[0], dna11[-1]
						startA2, endA2 = dna12[0], dna12[-1]
						startB1, endB1 = dna21[0], dna21[-1]
						startB2, endB2 = dna22[0], dna22[-1]

						description_string = "<p>{0}.{1} -> A</p>\n<p>{0}.{2} -> B</p>\n<p>{3}.{4} -> C</p>\n<p>{3}.{5} -> D</p>".format(code1, dna_chainA1, dna_chainA2, code2, dna_chainB1, dna_chainB2)
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
						
						rename("../tmp/StructAlign/{2}/{0}@_{1}@.pdb".format(code1, code2, random_name), "../tmp/StructAlign/{2}/{0}{3}_{1}{4}.pdb".format(code1, code2, random_name, chain1, chain2))
						
						print ("<center><H2>Your job was completed successfully</H2></center>\n")
						print '<table border="1" width="100%">\n<TR>\n<td>\n'
						print ( "<p>Structures {} chain {} and {} chain {} were aligned with score: {}</p>".format(code1, chain1, code2, chain2, score) )
					
						print ( '<p>To download the structural alignment follow this <a href="../tmp/StructAlign/{}/{}{}_{}{}.pdb">link</a></p> '.format(random_name, code1, chain1, code2, chain2) )


						#system("chmod a=rw ../tmp/StructAlign/{}{}_{}{}.pdb".format(code1, chain1, code2, chain2))
						print "<p>The nucleotides with max measure: {}.{}:{}, {}.{}:{}</p>".format(code1, dna_chainA1, maxA, code2, dna_chainB1, maxB)
						print ("<p>The alignment of DNAs:</p>" )
						middle1 = ''
						for index, nt in enumerate(dna1_chain1):
							if dna2_chain1[index] == nt:
								middle1 += '|'
							else:
								middle1 += ' '
						middle2 = ''
						for index, nt in enumerate(dna1_chain2):
							if dna2_chain2[index] == nt:
								middle2 += '|'
							else:
								middle2 += ' '

						print "<pre>{}.{}[{:>4}-{:<4}]: {}     {}.{}[{:>4}-{:<4}]: {}\n\t\t   {}     {:19}{}\n{}.{}[{:>4}-{:<4}]: {}     {}.{}[{:>4}-{:<4}]: {}</pre>".format(code1, dna_chainA1, startA1, endA1, dna1_chain1, code1, dna_chainA2, startA2, endA2, dna1_chain2, middle1, ' ', middle2, code2, dna_chainB1, startB1, endB1, dna2_chain1, code2, dna_chainB2, startB2, endB2, dna2_chain2)  
						#print "<pre>{}.{}[{:>4}-{:<4}]: {}\n\t\t   {}\n{}.{}[{:>4}-{:<4}]: {}</pre>".format(code1, dna_chainA2, startA2, endA2, dna1_chain2, middle2, code2, dna_chainB2, startB2, endB2, dna2_chain2)  

						print '      \n<script type="text/javascript"> \n \
		jmolSetCallback ( "messageCallback", "jmMessage" ); \n\
		jmolApplet ([600, 600], "load ../tmp/StructAlign/{}/{}{}_{}{}.pdb; cartoon only; color chain; background white" ); \n\
			jmolBr(); \n\
			jmolMenu(["background white","background black"]); \n\
	      </script> \n'.format(random_name, code1, chain1, code2, chain2)
					
						print '</td>\n<td width="50%" style="vertical-align: top">\n'
						print "<center><H3>Description</H3></center>"
						print warnings
						print '<p>To avoid possible overlapping of chain names, they were changed. The resulting pdb-file contains new names!</p>'
						print "<p>DNA chains:</p>"
						print description_string
						"""print "<p>{}.{} -> A</p>".format(code1, dna_chainA1)
						print "<p>{}.{} -> B</p>".format(code1, dna_chainA2)
						print "<p>{}.{} -> C</p>".format(code2, dna_chainB1)
						print "<p>{}.{} -> D</p>".format(code2, dna_chainB2)"""
						print "<p>Protein chains:</p>"
						print "<p>{}.{} -> E</p>".format(code1, chain1)
						print "<p>{}.{} -> F</p>".format(code2, chain2)
						print "</td></TR>"
						break
					
  
			else:
				print ( "<H3>Error</H3>" )
				print ( "<p>PDB entry <B>"+code2+"</B> does not exist</p>" )

print "</BODY>\n</HTML>"
		
system("/home/popov/bin/clean_temp.sh")

