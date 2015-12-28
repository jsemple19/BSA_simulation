from recombSim2 import *
import numpy as np
import os
import time
import random

mySeed=2015
random.seed(mySeed)

loci=15000
mutRate=3.1e-09
#chrRegions=np.array([0, 0.5, 4.5, 12.5, 16.5, 17]*loci/17,dtype=float)
#regRecRates=np.array([0, 7, 1, 7, 0]/15, dtype=float)

#qtls1=np.array([0.2*loci],dtype=int)
#qtls2=np.array([0.5*loci],dtype=int)
#qtls3=np.array([0.15*loci, 0.5*loci],dtype=int)
#qtls4=np.array([0.05*loci, 0.2*loci],dtype=int)
#qtls5=np.array([0.35*loci, 0.5*loci],dtype=int)
#qtls6=np.array([0.2*loci, 0.35*loci],dtype=int)


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
#selectDepths=[1,5,10,20,30,40,50]
#seqDepths=[25,50,100,150,200,250]
#qtls=[qtls1,qtls2,qtls3]
#effectSize=[1]
useRLE=True


baseDir=os.getcwd()+'/pySim_'+time.strftime("%Y%m%d")+'_'+str(loci)+'loci_'+str(numSims)+'sims/'
try:
    os.stat(baseDir)
except:
    os.mkdir(baseDir)       


attributes='q{}_N{}_G{}_sl{}_sq{}.txt'



for N in Ns:
    for numSim in range(0,numSims):
        population=Population(N)
        population.createF2(founder1,founder2,mutRate,RR)
        for G in Gs:
            population.createRIpop(G)
            writePopToFile(population,baseDir,numSim)

#fig=plt.plot(range(loci),alleleFreq1,'r',range(loci),alleleFreq,'b')
#plt.axis([0,loci,0,1])
#plt.xlabel('loci')
#plt.ylabel('allele frequency')
#plt.suptitle('whole genome vs per-locus selection for seq',fontsize=13)
#plt.vlines(qtls3,0,1,colors='0.4',linestyles='dashed')
#plt.hlines(0.5,0,loci,colors='0.4',linestyles='dotted')
#plt.title('N=1000, gen=5, selectDepth=0.1, seqDepth=100',{'fontsize':9})
#plt.figlegend(fig,('per-locus','whole-genome'),loc='lower right', 
#              bbox_to_anchor=(0.85, 0.2),fontsize=10)
#plt.savefig('genomeVslocus.pdf',papertype='a4',format='pdf')