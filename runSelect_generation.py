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

qtls1=np.array([0.2*loci],dtype=int)
qtls2=np.array([0.5*loci],dtype=int)
qtls3=np.array([0.15*loci, 0.5*loci],dtype=int)
#qtls4=np.array([0.05*loci, 0.2*loci],dtype=int)
#qtls5=np.array([0.35*loci, 0.5*loci],dtype=int)
#qtls6=np.array([0.2*loci, 0.35*loci],dtype=int)


#create vector of recombination rate along chr specific for C. elegans
recRate=celegRR(loci)
RR=recRate.perLocusRecRates()

# create genome template
genome3=haploidGenome(numLoci=loci, mutRate=mutRate, recProb=RR, useRLE=False)

#founder1=founderGenome(genome3)
#founder2=founderGenome(genome3)

#founder1.createFounder(genotype=0)
#founder2.createFounder(genotype=1)

numSim=50#sys.argv[1]
Ns=[10000]#[100,500,1000,5000,10000,50000,100000]
Gs=[5,10,15,20,25,30]#[5,10,15,20,25,30]
selectDepths=[10]#[1,5,10,20,30,40,50]#[1,5,10,20,30,40,50]
seqDepths=[100]#[25,50,75,100,125,150,200,250]
qtls=[qtls1,qtls2,qtls3]
effectSize=[1]


baseDir=os.getcwd()
inDir=baseDir+'/populations_pysim_'+str(loci)+'loci_'+str(numSim)+'sims_Dec2015/'
outDir=baseDir+'/pySim'+time.strftime("%Y%m%d")+'_'+str(loci)+'loci_'+'TgenNum/'
try:
    os.stat(outDir)
except:
    os.mkdir(outDir) 



attributes='q{}_N{}_G{}_sl{}_sq{}_sm{}.txt'


for N in Ns:
    for G in Gs:
        for sim in range(numSim):
            filename=inDir+'pop_N{}_G{}/sim{}.npz'.format(N,G,sim)
            population=readPopFromFile(filename,genome3) 
            for selectDepth in selectDepths:
                control=population.randomSelect(selectDepth) 
                for q in range(len(qtls)):                       
                    selected=population.pheSelect(qtls[q], effectSize, selectDepth)
                    for seqDepth in seqDepths:
                        contFreq=control.seqSample(seqDepth)
                        selectFreq=selected.seqSample(seqDepth)
                        values=[q,N,G,selectDepth,seqDepth,numSim]
                        writeToTxtFile(contFreq,selectFreq,outDir,
                                       attributes,values)
                                       
with open(outDir+'logfile.txt','at') as f:
    f.write(time.strftime("%Y%m%d")+'\n')
    f.write("number of simulations="+str(numSim)+'\n')
    f.write("size of populations="+' '.join(['%g' %i for i in Ns])+'\n')
    f.write("number of generations="+' '.join(['%g' %i for i in Gs])+'\n')
    f.write("depth of selection(%)="+' '.join(['%g' %i for i in selectDepths])+'\n')
    f.write("sequencing depth ="+' '.join(['%g' %i for i in seqDepths])+'\n')
    f.write("qtl positions=")
    for q in qtls:
        f.write(' '.join(['%g' %i for i in list(q)])+', ')
    f.write("\n")
    f.write("effect size="+' '.join(['%g' %i for i in effectSize])+'\n')
    f.write(__file__+'\n')
    f.write("set seed to="+str(mySeed)+'\n')