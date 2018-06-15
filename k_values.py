from odbAccess import *
from time import sleep
import math
from numpy import *
from sys import argv
from grow import getAllNodes, bonded, Cs, L0, Ltot


steps = (Ltot-L0)/2+1
thickness =[float(i) for i in argv[1:-1]]
name = argv[-1]


def getF( jobName):
	""" Extract the reaction force at crack tip 
	from the result file."""
	
	name = jobName + '.odb'
	odb = openOdb(name)
	F = fabs(odb.steps['Load'].frames[-1].fieldOutputs['RF']
		.values[0].data[1])
	return F
	
	
def getU(jobName,formertip):
	""" Extract node displacement after the crack tip
	from the result file."""
	
	name = jobName + '.odb'
	odb = openOdb(name)
	for node in odb.steps['Load'].frames[-1]\
	            .fieldOutputs['U'].values:
		if node.nodeLabel==formertip:
			U =  fabs(node.data[1])
	return U
	
	
def G(NodeF, NodeDisp, L):
	""" Calculate the energy releasing rate."""
	
	global Cs, thickness	
	# determine the local thickness at the crack tip:
	i = 0
	for C in Cs:
		if L < C:
			break
		else:
			i += 1
	th = list(thickness)
	Ns = len(thickness)
	for n in range(Ns):
		th.append(thickness[Ns-1-n])
	t = th[i]
	G = NodeF*NodeDisp/t
	return G
	
	
def K(G):
	""" Calculate the stress intensity factor."""
	
	K = (G*72400)**0.5
	return K
	

allnode = getAllNodes
bondedset = bonded(L0, allnode)

K_values = open('K_values_%s.txt' % name, 'w')
count = 0
for L in range(L0, Ltot+1):
	jobName = "%s_Crack%d" % (name, L)
	formertip = bondedset.pop(0)
	if count == 0 and L%2 ==0:
		for C in Cs:
			if C - L == 1:
				count = 2
			if C - L == 2:
				count = 4
	if count or L%2 == 0:
		if count != 0:
			count -= 1
		NodeDisp = getU(jobName, formertip)

		print >> K_values, "%10.6f" % L,
		NodeF    = getF(jobName)
		print >> K_values, "%12.6f" % NodeF,
		print >> K_values, "%12.6f" % NodeDisp,
		G_VCCT = G(NodeF, NodeDisp,L)
		K_VCCT = K(G_VCCT)
		print >> K_values, "%12.6f" % G_VCCT,
		print >> K_values, "%12.6f" % K_VCCT

K_values.close()
