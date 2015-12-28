from recombSim2 import *
import numpy as np
import os
import time
import random

mySeed=2015
random.seed(mySeed)

loci=15000
mutRate=3.1e-09

#create vector of recombination rate along chr specific for C. elegans
recRate=celegRR(loci)
RR=recRate.perLocusRecRates()

# create genome template
genome3=haploidGenome(numLoci=loci, mutRate=mutRate, recProb=RR, useRLE=True)

founder1=founderGenome(genome3)
founder2=founderGenome(genome3)

founder1.createFounder(genotype=0)
founder2.createFounder(genotype=1)

numSims=10
Ns=[5000]#[100,500,1000,5000,10000,50000,100000]
Gs=[5,10,15,20,25,30]
useRLE=True


baseDir=os.getcwd()+'/pySim_'+time.strftime("%Y%m%d")+'_'+str(loci)+'loci_'+str(numSims)+'sims/'
try:
    os.stat(baseDir)
except:
    os.mkdir(baseDir)       

for N in Ns:
    for numSim in range(0,numSims):
        population=Population(N)
        population.createF2(founder1,founder2,mutRate,RR)
        for G in Gs:
            population.createRIpop(G)
            writePopToFile(population,baseDir,numSim)

