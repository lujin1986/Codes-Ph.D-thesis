#!/usr/bin/python
from pylab import *
from deap import base, creator, tools, algorithms
from scipy import interpolate
import pickle
from os import system, access
from time import sleep
import glob
from Queue import Queue, Empty
from threading import Thread
from random import randint

indList = Queue() # set up the queue for parallization of jobs
Nind=20           # number of individuals in the population
NGEN=40           # maximum number of generations
CXPB = 0.5        # crossover rate
MUTPB = 0.2       # mutation rate at individual level
Glength = 19      # number of digits in the genotype of an individual
digits = 3        # number of digits representing the thickness of a 
				  # section
ths = 1           # number of threads
s = Glength

# see detailed explanation of base, creator and toolbox method 
# in DEAP documentation:
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
toolbox = base.Toolbox()
# attribute generator:
toolbox.register("attr_bool", randint, 0, 1)
# randomly generate the genotype:
toolbox.register("get_indi", tools.initRepeat, creator.Individual, 
    toolbox.attr_bool, Glength)

# setting up important GA parameters:
toolbox.register("mate", tools.cxOnePoint)                       
toolbox.register("mutate", tools.mutFlipBit, indpb=0.2)
toolbox.register("select", tools.selTournament, tournsize=3)	

def clearup():
	"""Delete the files that are generated from simulation."""
	
	Ext = ['prt', 'msg', 'sta', 'com', 'sim', 'odb', 
			'inp', 'dat', 'log']
	for ext in Ext:
		cmd = "rm -f *.%s" % ext 
		system(cmd)
		
		
def individual_():
	""" The method to randomly generate a genotype.""" 
	
	individual = toolbox.get_indi()
	return individual

	
def population_(n):
	"""The method to generate a randomized population with no duplication."""
	
	pop = []
	for i in range(n):
		if i == 0:
			pop.append(individual_())
		else:
			new = individual_()
			# if "new" has appeared, regenerate "new":
			flag = 0
			for ind in pop:
				if ind[:s] == new[:s]:
					flag = 1
			while flag:
				new = individual_()
				flag = 0
				for ind in pop:
					if ind[:s] == new[:s]:
						flag = 1
			pop.append(new)
	return pop	
	
# register the two methods in toolbox:	
toolbox.register("individual", individual_)
toolbox.register("population", population_)


def grow(List):
	"""Automated generation of FEM models
	and execuation of simulations."""
	
	# get the queue of jobs for parallization:
	global indList
	
	# execute jobs until the queue is empty:
	while True:
		try:
			job=indList.get(block=False)
		except Empty:
			break
		else:
		    # define the name in the form: 
			# "index of generation"_"index of individual":
			name = "%d_%d" % (job[1], job[2]) 
			thickness = decoder(job[0])
			flag2 = 0
			count = 0
            # check if this individual has been evaluated:
			in this generation:
			for ind in List[:job[2]]:
				if ind == job[0]:
					flag2=1
					break
				else: count+=1
				
			# if it has been evaluated, 
			# directly copy the extracted K profile:
			if flag2:
				while not os.access(sta_file_name, os.R_OK):
					sleep(1)
				sleep(10)
				cmd = "cp K_values_%d_%d.txt K_values_%s.txt" \
				% (job[1], count, name)
				system(cmd)
			else:
				# call the subroutine grow.py to 
				# generate FEM model and to run the simulation:
				print "start growing %s" % name
				cmd = "python grow.py"
				for t in thickness:
					cmd += ' %s' % str(t)
				cmd += ' %s' % name
				print "submitting job: %s" % cmd
				system(cmd)
		
		
def getN(gen,index,c_, m_):
	"""Calculated the fatigue life 
	based on the extracted K profile."""
	
	name = "%d_%d" % (gen, index)
	filename = 'K_values_' + name + '.txt'
	data = loadtxt(filename)
	a = data[:,0]
	K = data[:,-1]*0.001**0.5
	dadN = c_*(K)**m_
	dN = 1/dadN
	x_new = np.linspace(6, 145, 10000)
	tck = interpolate.splrep(a, dN)
	dN_new = interpolate.splev(x_new, tck)
	# integration dN/da over a from 6 mm to 145 mm:
	N = np.trapz(dN_new, x_new)
	return (N, )

	
def decoder(individual):
	"""Decode the genotype into a list of thicknesses
	of different sections."""
	
	global digits, scheme
	# decode binary genotype firstly into a set of integers:
	gene = individual[:-4]
	integers = []
	for i in range((len(individual)-4)/digits):
		j = 0
		inte = 0
		while j < digits:
			inte += gene[i* digits+j]*2**(digits-j-1)
			j += 1
		integers.append(inte)
		
	# decode the set of integers into the set of thickness values 
	# considering both the boundaries 
	# and constant weight constraints: 
	thickness = []
	for inte in integers:
		thickness.append(1.9+inte/7.0*2.25)
	M = 2.9*len(thickness)
	tot = sum(thickness) 
	# small compensation term "com" added to each thickness section
	# to maintain constant weight:
	com = (M - tot)/len(thickness)
	for t in range(len(thickness)):
		thickness[t] += com
	while True:
		flag = 0
		i=0
		for t in thickness:
			if  t < 1.899999:
				comp = 1.9-t
				for j in range(len(thickness)):
					if j != i:
						thickness[j] -= comp/(len(thickness)-1)
				thickness[i] = 1.9

				flag =1
			elif t > 4.150001:
				comp = t - 4.15
				for j in range(len(thickness)):
					if j != i:
						thickness[j] += comp/(len(thickness)-1)
				thickness[i] = 4.15
				flag =1
			i+=1
		if flag==0:		
			break
	return thickness
	
	
def history(pop, append = True):
	"""Record the present generation and its relevant statistics
	   in a plain text file."""
	   
	# gather all the fitnesses and phenotypes in one list:
	fits = [ind.fitness.values[0] for ind in pop]
	thicknesses = [decoder(ind) for ind in pop]
	final = [i+[j] for i, j in zip(thicknesses, fits)]		
	savetxt("life_Gen%d.txt" % initGEN, array(final))
	
	# produce the statistic report:
	length = len(pop)
	mean = sum(fits) / length
	sum2 = sum(x*x for x in fits)
	std = abs(sum2 / length - mean**2)**0.5
	ind = [i for i, j in enumerate(fits) if j==max(fits)]    
	fittest = pop[ind[0]]	
	if append:
		history = open("history.txt", 'a')
	else:
		history = open("history.txt", 'w')
	print >> history, "====== Generation %d results ======" \
	% initGEN
	print >> history, "  Min %s" % min(fits)
	print >> history, "  Max %s" % max(fits) 
	print >> history, "  Avg %s" % mean 
	print >> history, "  Std %s" % std 
	print >> history, "genotype of fittest individual: %s" \
	% fittest 
	print >> history, "phenotype of fittest individual: %s"\
	% decoder(fittest) 
	history.close()
	
	
def parallization(work):
	""" Parallilization of jobs."""
	
	workers = [Thread(target=grow, args=(work,)) 
	for i in range(ths)]
		for worker in workers:
			worker.start()
		for worker in workers:
			worker.join()	

	
def main(restart = False, elitism = True):
	""" The main routine to perform FEM-GA coupled optimization.
	    restart: to restart the optimization from a broken point,
		elitism: the switch on/off of the elitism in the optimization.
	"""
	
	initGEN = 0    # index of initial generation     
	if restart:
		clearup()
		# load the archive for evaluated individuals:
		f = open("valid_ind.txt" , 'rb')    
		valid_ind = pickle.load(f)
		f.close()
		# list of files containing all the genotypes 
		# in the beginning of each generation:
		nameListOff = glob.glob("offspring_Gen_*.txt")  
		# list of files containing all the genotypes and 
		# their fitnesses at the end of each generation:
		nameListPop = glob.glob("population_Gen_*.txt")
		a =  max([int(i[14:-4]) for i in nameListOff])
		b =  max([int(i[15:-4]) for i in nameListPop])
		if b < a:
			initGEN = a
			print "restart at generation %d " % initGEN
			f = open("offspring_Gen_%d.txt" % initGEN, 'rb')
			offspring = pickle.load(f)
			f.close()
			# only evaluate the individuals 
			# with an invalid fitness value:
			invalid_ind = [ind for ind in offspring 
						   if not ind.fitness.valid]
			
			# filter out the individuals that have been validated
			# according to the achive of all evaluated individuals:
			for indin in invalid_ind:
				for indv in valid_ind:
					if indin == indv:
						indin.fitness.values = indv.fitness.values
			invalid_ind = [ind for ind in invalid_ind 
						   if not ind.fitness.valid]	
			
	        # calculate the fitness values for the individuals, 
			# the simulations of which have been performed:
			Klist = glob.glob("K_values_%d_*.txt" % initGEN)
			finish_n = max([int(i[len("K_values_%d_"  % initGEN):-4])
							for i in Klist])
			for index in range(0, finish_n+1):
				invalid_ind[index].fitness.values= \
				getN(initGEN, index,  4e-7, 2.4)
				valid_ind.append(invalid_ind[index])
			
			# execute the FEM simulations 
			# for the rest of population in multiple threads:
			rest = range(finish_n+1, Nind)
			if rest:
				for index in rest:
					indList.put([invalid_ind[index], 
								 initGEN, index])
		    parallization(invalid_ind)
			for index in rest:				
				invalid_ind[index].fitness.values= \
				getN(initGEN, index,  4e-7, 2.4)
				valid_ind.append(invalid_ind[index])	
				
			# replace the old population with 
			# the newly evolved population:		
			pop[:]=offspring
			# pickle all the genotypes and 
			# the corresponding fitness values of the population:
			f = open("population_Gen_%d.txt" % initGEN, 'wb')
			pickle.dump(pop, f)
			f.close()
			history(pop)
			initGEN = initGEN+1
			clearup()
		else: 
			f = open("population_Gen_%d.txt" % b, 'rb')
			pop = pickle.load(f)
			f.close()
			initGEN += 1
			

	for g in range(initGEN,NGEN):	
		if g==0:
			print "starting generation %d."	% g
			#randomly generate the initial population:
			pop = toolbox.population(n=Nind)
			valid_ind = []
			# pickle the state of the population 
			# in the beginning of the generation:
			f = open("offspring_Gen_0.txt", 'wb')
			pickle.dump(pop, f)
			f.close()
			offspring = pop
			print "start Gen %d" %g
			for index in range(Nind):
				indList.put([offspring[index], g, index])
			parallilization(offspring)
			fitnesses = [getN(g, index, 4e-7, 2.4) 
						 for index in range(Nind)] 
			for ind, fit in zip(offspring, fitnesses):
				ind.fitness.values = fit
				valid_ind.append(ind)
			f = open("population_Gen_%d.txt" % g, 'wb')
			pickle.dump(offspring, f)
			f.close()
			pop[:] = offspring
			fittest = history(pop, append=false)
		
		else:
			print "Generation %d is being generated..." % g
			# perform selection based on fitness:
			offspring = toolbox.select(pop, len(pop))
			# clone the selected individuals:
			offspring = list(map(toolbox.clone, offspring))
			
			# apply crossover and mutation on the offspring:
			for child1, child2 in zip(offspring[::2], \
									  offspring[1::2]):
				if random() < CXPB:
					toolbox.mate(child1, child2)
					del child1.fitness.values
					del child2.fitness.values
			for mutant in offspring:
				if random() < MUTPB:
					toolbox.mutate(mutant)
					del mutant.fitness.values
					

			f = open("offspring_Gen_%d.txt" % g, 'wb')
			pickle.dump(offspring, f)
			f.close()
			if elitism:
				if fittest not in offspring:
					offspring[0] = fittest
			invalid_ind = [ind for ind in offspring 
						   if not ind.fitness.valid]
			for indin in invalid_ind:
				for indv in valid_ind:
					if indin == indv:
						indin.fitness.values = indv.fitness.values
			invalid_ind = [ind for ind in invalid_ind 
						   if not ind.fitness.valid]			
			for index, ind in enumerate(invalid_ind):
				indList.put([ind, g, index])
			parallilization(invalid_ind)
			fitnesses = [getN(g, index, 4e-7, 2.4) 
			             for index in range(len(invalid_ind))] 
			for ind, fit in zip(invalid_ind, fitnesses):
				ind.fitness.values = fit
				valid_ind.append(ind)
			pop[:] = offspring
			f = open("population_Gen_%d.txt" % g, 'wb')
			pickle.dump(pop, f)
			f.close()
			fittest = history(pop)
		clearup()
		
		# update the achive for evaluated individuals:
		f = open("valid_ind.txt" , 'wb')
		pickle.dump(valid_ind,f)
		f.close()
		
if __name__ == "__main__":
	main(restart=0)	
