from odbAccess import *
from time import sleep
import math
import numpy as np
from sys import argv
from Queue import Queue, Empty
from threading import Thread

# positions of borders between different thickness sections:
Cs =  [19,33,47,61,75,89,103,117,131]       
Nworker = 2			# number of threads
joblist = Queue()  	# set up the queue for parallization of jobs

# receive the thickness values passed from main.py:
thickness = [float(a) for a in argv[1:-1]]
name = argv[-1]
	
L0 = 6 #mm                                  initial crack length                       
Ltot = 145 #mm                              final crack length
steps = (Ltot-L0)/2+1 #                     number of increments

				
def getAllNodes():
	""" Get the list of all nodes in the model.
	    Each Node in the list has the format: 
		["node label", "x coordinate", "y coordinate"].
	"""
	
	f = open('other/model.inp', 'r')
	allnode = []
	flag = 0
	for line in f.readlines():
		data = line.split()
		if data[0] == '*Node':
			flag = 1
			continue
		if flag:
			if data[0][0] != '*':
				allnode.append([float(i) for i in line.split(',')])
			else: break	
	f.close()
	return allnode
		

def bonded(L0, Nlist):
	""" Get the label lists of nodes located on the crack line.
		L0 is the initial crack length.
		Nlist is the list of all nodes.
	"""
	tol = 0.01
	NodeBonded = [] # all the nodes along the crack line
	NodeLabelBonded= [] # label of all the nodes on the crack line	
	for N in Nlist:
		if fabs(N[2]) < tol: 
			if N[1] > L0-1 - tol:
				NodeBonded.append(N)
	# sort the node list according to the x coordinate:
	NodeBonded = NodeBonded.sort(key=lambda x: x[1]) 
	NodeLabelBonded = [int(i[0]) for i in NodeBonded]
	return NodeLabelBonded	
	
	
def ABQinp(bondedSet,pretip, jobName, name):
	"Prepare input file for simulation"
	
	symmynodes = ' '
	for i in bondedSet:
		if ind%8 == 0:
			symmynodes += str(i)+ ',\n'
		else:
			symmynodes += str(i)+ ', '
	f = open(jobName+'.inp', 'w')
	model = open(name+'.inp', 'r')
	lines = model.readlines()
	for line in lines:
		print >> f, line.replace('symmynodes', symmynodes)\
		.replace('pretipnode', str(pretip))\
		.replace('tipnode', str(bondedSet[0])),
	f.close()
	

def runSim(i):
	""" Run the simulation jobs that has been put into the queue. """
	
	global joblist
	while True:
		try:
			jobName=joblist.get(block=False)
		except Empty:
			break
		else:	
			if not access("%s_min.odb" % jobName, R_OK):
				cmd = 'abaqus job='+jobName+' interactive'
				system(cmd)


def growup(name):
	""" Run the series of simulation jobs with the crack 
		incrementally extending from an initial length to 
		a final length. And extract the stress intensity
		factor profile from the result files of the simulations.
	"""
	
	# prepare to write the common input file
	# for the series of simulation:
	f = open(name+'.inp', 'w')
	# read the content from the template of the input file:
	model = open('other/model.inp', 'r')
	lines = model.readlines()
	
	# write the thicknesses of the crenellation pattern 
	# into corresponding positions of input file:	
	markers = ['th%d' % n for n in range(1,len(thickness)+1)]
	for line in lines:
		for t, marker in zip(thickness, markers):
			line = line.replace(marker, str(t))
		print >> f, line,
	f.close()   

	allnode = getAllNodes()
	bondedset = bonded(L0, allnode)

	# the crack extends with an increment of 2 mm,
	# and if the crack tip is in the area near the crenellation steps
    # the increment is reduced to 1 mm
	# due to the sharp change of K values locally:
	count = 0
	for L in range(L0, Ltot+1):
		jobName = "%s_Crack%d" % (name, L)
		formertip = bondedset.pop(0)
		if count == 0 and L%2 == 0:
			for C in Cs:
				if C - L == 1:
					count = 2
				if C - L == 2:
					count = 4
		if count or L%2 == 0:
			ABQinp(bondedset, formertip, jobName, name)	
			joblist.put(jobName)
			if count != 0:
				count -= 1
	# parallization of simulations:
	workers = [Thread(target=runSim, args = (i,)) 
			   for i in range(Nworker)]
	for worker in workers:
		worker.start()
	for worker in workers:
		worker.join()

	# check if the last simulation job is finished:	
	sta_file_name =  '%s_Crack144.sta' % name
	while not access(sta_file_name, R_OK):
		sleep(1)
		
	# extracting K values from the result files 
	# by calling subroutine K_values.py: 
	print "start extracting Ks for %s" % name
	cmd = "abaqus  python getKs.py"
	for t in thickness:
		cmd += ' %s' % str(t)
	cmd += ' %s' % name
	system(cmd)		

	
growup(name)
