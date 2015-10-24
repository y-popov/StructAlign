#!/usr/bin python
current = "All"

shifts = []

g = open("random.txt", 'r')
measures1 = []
for line in g:
	measures1.append(float( line.split()[3] ))
measures1.sort()
g.close()

f = open("{}.txt".format(current), 'r')
measures2 = []
for line in f:
	measures2.append( float(line.split()[3]) )
measures2.sort()
f.close()

assay = []
m_assay = []
i=j=0
for k in range(len(measures1)+len(measures2)):
	if (measures1[i] < measures2[j]):
		assay.append(measures1[i])
		m_assay.append(1)
		i += 1
	else:
		assay.append(measures2[j])
		m_assay.append(2)
		j += 1
	if i==len(measures1):
		assay.extend(measures2[j:])
		m_assay.extend([2 for x in range(len(measures2)-j)])
		break
	elif j==len(measures2):
		assay.extend(measures1[i:])
		m_assay.extend([1 for x in range(len(measures1)-i)])
		break

def BSearch(assay, m_assay):
    return BSearch_partial(0, len(assay), assay, m_assay)

def BSearch_partial(start, end, array, m_array):
    right = end
    left = start
    i = start
    while True:
	i = (left+right)/2
	print 'I am in', array[i], 'position=', i
	lower_i = 0
	greater_i = 0
	for j in xrange(len(array)):
		if (array[j] < array[i]) and (m_array[j] == 2):
			lower_i += 1
		if (array[j] > array[i]) and (m_array[j] == 1):
			greater_i += 1
	print "Lower i: {}, Greater i: {}\n".format(lower_i, greater_i)
	if lower_i == greater_i:
	    return array[i]
	if lower_i > greater_i:
	    right = i
	elif lower_i < greater_i:
	    left = i

meow = BSearch(assay, m_assay)
shifts.append(meow)

print shifts

