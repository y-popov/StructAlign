import re

list = []
d = {}

f = open("Homeobox_output.txt", 'r')
for line in f:
	temp = re.findall(r'new\/(.+).pdb-.+\/(.+).pdb\t(.+)', line)[0]
	if temp[0] not in d:
		d[temp[0]] = {}
		list.append(temp[0])
	d[temp[0]][temp[1]] = temp[2].split()[-1]


f.close()

#d['A_1jgg']['B_1du0'] = 'seg_f'

print ' '*6,
for i in list:
	print "{:6} ".format(i),
print ""

for i in range(len(list)):
	print list[i],
	for j in range(len(list)):
		if j>=i:
			s = d[list[i]][list[j]]
			if s.replace('.', '').isdigit():
				print "{:>6.0f} ".format(float(s)),
			else:
				print "{:>6} ".format(s),
		else:
			print " "*7,
	print ""

