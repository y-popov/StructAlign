#!/usr/bin/python

from sys import stdin
from os import system, access, F_OK, rename, makedirs
from random import choice
from string import ascii_uppercase, digits
from cgi import FieldStorage


def complement(x):
	alph = 'ATGC-?'
	compl_alph = 'TACG--'
	s = ''
	for i in x:
		s += compl_alph[alph.index(i)]
	return s


def StructAlign(pdb1, pdb2, outfile, warns, pairs, maxMs, ranges):
	code1 = pdb1[pdb1.rfind('/')+1:-1]
	code2 = pdb2[pdb2.rfind('/')+1:-1]
	
	chain1 = chain1_old = pdb1[-1].upper()
	chain2 = chain2_old = pdb2[-1].upper()
	
	range1 = ranges[code1][chain1] 
	range2 = ranges[code2][chain2] 
	
	pdb1_name = "../tmp/"+code1+'.pdb'
	pdb2_name = "../tmp/"+code2+'.pdb'
	
	output = code1+'@_'+code2+'@'
	
	try:
		system('StructAlign/align {0} {1} ../tmp/StructAlign/{2}/All/{3}.pdb {4} {5} {6} {7} {8} {9} /var/www/tools/tmp/StructAlign/{2}/result.txt 1 >> ../tmp/StructAlign/{2}/log.txt'.format(pdb1_name, pdb2_name, outfile, output, chain1, chain2, range1[0], range1[1], range2[0], range2[1]))
	except Exception as e:
		error.append(e)
	
	
	max_score = open("/var/www/tools/tmp/StructAlign/{}/result.txt".format(outfile), 'r')
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

			pdb_name = '../tmp/StructAlign/'+outfile+'/All/'+pair+'.pdb'
			rename('../tmp/StructAlign/'+outfile+'/All/'+output+'.pdb', pdb_name)
			
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

def make_latex_matrix(d, list):
	s = "\\(\\begin{matrix}"
	leng = len(list[0])

	for i in range(len(list)):
		s += "\\textbf{"+list[i]+"} & "
		for j in range(len(list)):
			if i>=j:
				n = d[list[i]][list[j]]
				s += "{:>7.0f} &".format(float(n))
		s += "\\\\"
	s += ' & '
	for i in list:
		s += "\\textbf{"+i+"} & "
	s += "\\end{matrix}\\)"
	return s

def make_table_matrix(d, list):
	s = "<table class='table table-hover table-condensed'> \n <tr>"
	leng = len(list[0])

	for i in range(len(list)):
		s += "<td><b>"+list[i]+"<b></td>"
		for j in range(len(list)):
			if i>=j:
				n = d[list[i]][list[j]]
				s += "<td>{:>7.0f}</td>".format(float(n))
			else:
				s += " <td></td>"
		s += "</tr>"
	s += ' <tr><td></td> '
	for i in list:
		s += "<td><b>"+i+"</b></td> "
	s += "</tr></table>"
	return s

def find_repr(m):
	mins = [min(x) for x in m]
	best = max(mins)
	i = mins.index(best)
	
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

	system("StructAlign/realign ../tmp/{}.pdb {} {} {} {} ../tmp/StructAlign/{}/All/{}.pdb {} {} {} {} {} >> ../tmp/StructAlign/{}/log.txt".format(repr_pdb[:-1], repr_pdb[-1], maxM1[0][1], maxM1[0][4], maxM1[1][1], random_name, pair, index, pdb_name, target[-1], chain1, chain2, random_name))
	
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

def print_fasta(fastas):
	for handle in fastas:
		print "{}.{}[{:>4}-{:<4}]: {}".format(handle.pdb, handle.chain, handle.num_seq[0], handle.num_seq[-1], '-'*handle.gaps[0]+handle.seq+'-'*handle.gaps[1])



print "Content-type: text/html\n\n"
print "<HTML>\n<HEAD>"
print ' <script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script> \
 <script type="text/javascript" src="/tools/jsmol/Jmol.js"></script> \n \
  <script type="text/javascript"> \n \
    jmolInitialize("/tools/jsmol/"); \n \
  </script> \n \
  ' 
print '''
	<meta charset="utf-8">
	<meta http-equiv="X-UA-Compatible" content="IE=edge">
	<meta name="viewport" content="width=device-width, initial-scale=1">
	<meta name="description" content="A program for alignment of structures of DNA-protein complexes">
	<meta name="author" content="Yaroslav Popov <syav.popoff@yandex.ru>">
	<link rel="shortcut icon" href="StructAlign/favicon.ico">
<TITLE>Result</TITLE>
	<link href='https://fonts.googleapis.com/css?family=Roboto:400,300,500,100' rel='stylesheet' type='text/css'>
	<link href="../Bootstrap3/css/bootstrap.min.css" rel="stylesheet">
	<link href="../StructAlign/StructAlign.css" rel="stylesheet" type="text/css">
	
</HEAD>\n

<style>
	.main {
		margin: 0 10px;
	}
</style>'''
print "<BODY>"
print '\
  <script type="text/javascript">    \n \
    function jmMessage ( apName, Value )  \n \
    {   \n \
      var el = document.getElementById ( "jmConsole" ); \n \
      el.innerHTML = el.innerHTML + "<br>" + Value;   \n \
    }   \n \
  </script> \n \
  '


print '''
	    <div class="navbar navbar-inverse navbar-static-top" role="navigation">
	      <div class="container">
		<div class="navbar-header">
		  <a class="navbar-brand" href="../StructAlign.html">StructAlign</a>
		</div>
		<div class="collapse navbar-collapse">
		  <ul class="nav navbar-nav navbar-right">
		    <li><a href="../StructAlign.html">Home</a></li>
		    <li><a href="../StructAlign/StructAlign_feedback.html">Feedback</a></li>
		    <li><a href="../StructAlign/StructAlign_help.html">Help</a></li>
		    <li><a href="../StructAlign/StructAlign_about.html">About</a></li>
		  </ul>
		</div><!--/.nav-collapse -->
	      </div>
	    </div>
	    
	    <div class="main">
'''

is_multi = False
c = []
form = FieldStorage()
fileitem = form['file']
if fileitem.filename:
	fileitem = fileitem.file.read().splitlines()
	if len(fileitem) > 2:
		is_multi = True
	else:
		print '''<div class="alert alert-danger"><strong>Error!</strong> There are less then 2 PDBs in your request file</div>'''
else:
	fileitem = None
	for item in form.__iter__():
		c.append(str(item)+'='+str(form.getvalue(item)))

c.sort()
system('echo "{}" > ../tmp/StructAlign/log.txt'.format(c))
random_name = ''.join(choice(ascii_uppercase + digits) for i in range(10))

for d in c:
	(name,value) = d.split("=")
	if name == "a":
		is_multi = True if value=="1" else False
	
	if not is_multi:
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


				
				makedirs("../tmp/StructAlign/{}".format(random_name))
				system("touch ../tmp/StructAlign/{}/result.txt".format(random_name) )
				system("chmod a=rwx -R ../tmp/StructAlign/{}".format(random_name))
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

if is_multi:
	from Bio import Phylo 
	from Bio.Phylo import TreeConstruction
	from subprocess import Popen, PIPE
	
	pdbs = []
	chains = []
	starts = []
	ends = []
	ranges = {}
	print "<h1>This option is under construction</h1><br>"
	
	if fileitem:
		for query in fileitem:
			pdb, buff = query.split(".")
			chain, buff = buff.split(":")
			start, end = buff.split("--")
			
			pdbs.append(pdb.lower())
			chains.append(chain.upper() if chain else '@')
			starts.append(start if start else "zero")
			ends.append(end if end else "inf")
			
	else:
		for i in c:
			(name,value) = i.split("=")
			if name.startswith("chain"):
				chains.append(value.upper() if value else '@')
			elif name.startswith("start"):
				starts.append(value if value else "zero")
			elif name.startswith("end"):
				ends.append(value if value else "inf")
			elif name.startswith("zpdb"):
				pdbs.append(value.lower())

	
	for index, code in enumerate(pdbs):
		filename = "/mnt/databanks/pdb/"+code[1:3].lower()+"/pdb"+code+".ent.gz"
		if access(filename, F_OK): # this PDB entry exists
			if not access("../tmp/"+code+".pdb", F_OK):
				system("gunzip -c " + filename + "> ../tmp/"+code+".pdb")
		else:
			print ( "<H3>Error</H3>" )
			print ( "<p>PDB entry <B>'"+code+"'</B> does not exist</p>" )
			exit(1)
		chain = chains[index]
		pdbs[index] = code + chain
		r = [starts[index], ends[index]]
		if code not in ranges:
			ranges[code] = {chain:r}
		else:
			ranges[code][chain] = r
	
	makedirs("../tmp/StructAlign/{}".format(random_name))
	makedirs('../tmp/StructAlign/{}/All'.format(random_name))
	system("touch ../tmp/StructAlign/{}/result.txt".format(random_name) )
	system("chmod a=rwx -R ../tmp/StructAlign/{}".format(random_name))
	
	warnings = ''
	errors = []
	pairs = []
	maxMs = {}
	PhyloM = []
	scores = {}
	scoresM = []
	s = 0
	leng = len(pdbs)
	total = (leng**2+leng)/2.0
	
	
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
			if score is not None:
				pdbs[i] = pdbs[i][:-1]+chain1
				pdbs[j] = pdbs[j][:-1]+chain2
				if pdbs[i] not in scores:
					scores[pdbs[i]] = {}
					PhyloM.append([])
					scoresM.append([])
				scores[pdbs[i]][pdbs[j]] = score
				scoresM[i].append(score)
			else:
				if pdbs[i] not in scores:
					scores[pdbs[i]] = {}
					PhyloM.append([])
					scoresM.append([])
				scores[pdbs[i]][pdbs[j]] = 0
				scoresM[i].append(0)
	
	scoresM = TreeConstruction._Matrix([x[x.rfind('/')+1:] for x in pdbs], scoresM)
	print "<pre>"
	print_matrix("Scores", scores, pdbs)
	print "</pre>"
	print "<h4>Scores:</h4>"
	print make_latex_matrix(scores, pdbs)
	print make_table_matrix(scores, pdbs)
	
	distances = {}
	for i in range(leng):
		distances[pdbs[i]] = {}
		for j in range(i+1):
			distances[pdbs[i]][pdbs[j]] = (scores[pdbs[i]][pdbs[i]]+scores[pdbs[j]][pdbs[j]])/2.0 - scores[pdbs[i]][pdbs[j]]
			PhyloM[i].append(distances[pdbs[i]][pdbs[j]])
	PhyloM = TreeConstruction._DistanceMatrix([x[x.rfind('/')+1:] for x in pdbs], PhyloM)
	print "<pre>"
	print_matrix("Distances", distances, pdbs)
	print "</pre>"
	print "<h4>Distances:</h4>"
	print make_latex_matrix(distances, pdbs)
	print make_table_matrix(distances, pdbs)
	
	tree = TreeConstruction.DistanceTreeConstructor().upgma(PhyloM)
	multi_tree = "../tmp/StructAlign/"+random_name+"/multi.tre"
	multi_ps = "../tmp/StructAlign/"+random_name+"/multi.ps"
	multi_png = "../tmp/StructAlign/"+random_name+"/multi.png"
	Phylo.write(tree, multi_tree, "newick")
	arguments = ["fdrawgram", multi_tree, multi_ps, "-previewer", "n"]
	fdrawgram = Popen(arguments, stdout = PIPE, stderr = PIPE)
	(fdgtrash, fdgerror) = fdrawgram.communicate()
	if fdrawgram.returncode == 0:
		arguments = ["gs", "-sDEVICE=pngmono", "-sOutputFile=" + multi_png, "-r300", "-dBATCH", "-dNOPAUSE", multi_ps]
		gs = Popen(arguments, stdout = PIPE, stderr = PIPE)
		(gstrash, gserror) = gs.communicate()
		if gs.returncode == 0:
			pass
	
	print "<pre>"
	Phylo.draw_ascii(tree)
	print "</pre>"
	print "<img src='../tmp/StructAlign/"+random_name+"/multi.png' width=500>"
	
	
	
	repres, repr_index = find_repr(scoresM)
	
	multy_pdb = "../tmp/StructAlign/"+random_name+"/multi.pdb"
	f = open(multy_pdb, 'w')
	f.write( "HEADER{:>60}\n".format("MultiStructAlign") )
	f.write( "TITLE{:>60}\n".format("title") )
	f.close()
	
	c = 1
	description = "\n"
	for pair in pairs:
		if repres in pair:
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
	
	multy_fasta1 = "../tmp/StructAlign/"+random_name+"/multi1.fasta"
	multy_fasta2 = "../tmp/StructAlign/"+random_name+"/multi2.fasta"
	g = open(multy_fasta1 ,'w')
	alignFASTA(fasta_align1, ref_fasta1)
	writeFASTA(g, fasta_align1)
	g.close()

	g = open(multy_fasta2 ,'w')
	alignFASTA(fasta_align2, ref_fasta2)
	writeFASTA(g, fasta_align2)
	g.close()
	
	multy_descr = "../tmp/StructAlign/"+random_name+"/multi.txt"
	f = open(multy_descr, 'w')
	f.write("Model	PDB_ID	Protein_Chain	DNA_Chain_1	DNA_Chain_2\n")
	for index, i in enumerate(description.split('\n')[1:-1]):
		buff = i.split()
		f.write("{}\t{}\t{}\t{}\t{}\n".format(buff[1][:-1], buff[2][:-1].upper(), buff[2][-1], fasta_align1[index].chain, fasta_align2[index].chain))
	f.close()
	
	
	print "<pre>"
	print_fasta(fasta_align1)
	print "</pre>"
	
	print "<pre>"
	print_fasta(fasta_align2)
	print "</pre>"
	
	system("zip -j ../tmp/StructAlign/{0}/multi ../tmp/StructAlign/{0}/multi.pdb ../tmp/StructAlign/{0}/multi.txt ../tmp/StructAlign/{0}/multi1.fasta ../tmp/StructAlign/{0}/multi2.fasta".format(random_name))
	
	print "<a href='../tmp/StructAlign/{}/multi.pdb'>PDB</a>".format(random_name)
	print "<a href='../tmp/StructAlign/{}/multi.zip'>ZIP</a>".format(random_name)
	
	
	print '      \n<script type="text/javascript"> \n \
		jmolSetCallback ( "messageCallback", "jmMessage" ); \n\
		jmolApplet ([600, 600], "load ../tmp/StructAlign/{}/multi.pdb; model 0; cartoon only; color chain; background white" ); \n\
			jmolBr(); \n\
	      </script> \n'.format(random_name)

	
		
print "</div></BODY>\n</HTML>"
		
system("/home/popov/bin/clean_temp.sh")

