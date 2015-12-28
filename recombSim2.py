# last modified 20151215 changed sequencing selection to use sampling with 
# replacement for all populations sizes

import numpy as np
import random as rnd
from itertools import groupby
import os
import random

#mySeed=2015
#random.seed(mySeed)


def encode(input_nparray):
    '''
    converts a numpy boolean array to numpy array with run length encoding
    array( [1,1,1,0,0,0], dtype=bool ) -> array( [[3,1],[3,0]], dtype=int )
    Note: using a list of lists might be faster but takes a bit more space
    not sure what the tradeoff is.
    '''
    return np.array([[len(list(rep)), int(val)] for val,rep in 
                    groupby(input_nparray)],dtype=int)

def decode(rle_lst):
    '''
    expands a numpy two dimensional array of run length encoding into a one
    dimensional numpy boolean array
    array( [[3,1],[3,0]], dtype=int ) -> array( [1,1,1,0,0,0], dtype=bool )
    '''
    return np.array(np.hstack([val] * rep for rep,val in rle_lst),dtype=bool)


class celegRR:
    '''
    object with C. elegans specific recombination rate data from Rockman et al.
    self.fractions takes the number of cM and transforms them to fraction of the
    total chromosome length
    self.rates takes the recombination rates in each of these regions, and adds 
    a small adjustment to that the tips of the chromosome have at least some low 
    recombination rate(>100x less then the centre of the choromosome)
    '''
    def __init__(self,numLoci=100):
        self.numLoci=numLoci
        self.fractions=np.array([[3.5, 22.1, 47.7, 25.4, 1.3],
            [2.0, 29.9, 46.7, 16.9, 4.5],
            [3.6, 23.4, 48.0, 20.9, 4.1],
            [4.1, 18.2, 51.9, 21.4, 4.4],
            [3.1, 25.1, 50.9, 18.1, 2.8],
            [3.2, 31.4, 35.8, 22.2, 7.4]],dtype=float)/100
        self.chrNames=["chrI","chrII","chrIII","chrIV","chrV","chrX"]
        self.adjust=np.vstack((np.array([0.01]*6),np.array([-0.01]*6),
                                        np.array([0]*6),np.array([-0.01]*6),
                                        np.array([0.01]*6))).transpose()
        self.rates=np.array([[0, 3.43, 1.34, 6.78, 0],
             [0, 4.92, 1.33, 8.47, 0],
             [0, 7.83, 1.17, 7.24, 0],
             [0, 7.65, 1.05, 3.64, 0],
             [0, 3.22, 1.32, 5.47, 0],
             [0, 3.81, 1.70, 5.14, 0]],dtype=float)+self.adjust


    def perLocusRecRates(self,chr="chrIII"):
        '''
        perLocusRecRates(str) -> nd.array
        Returns a numpy array of length numLoci with specific recombination 
        probability at each locus, such that they sum to a total probability is 
        1 over the whole chromosome (array). If none of the C. elegans 
        chromosomes are specified then an array of uniform probabilities
        will be returned.
        '''
        if (chr not in self.chrNames):
            probRec=np.ones(self.numLoci,dtype=float)/self.numLoci
            return probRec
        i=self.chrNames.index(chr)
        probRec=np.hstack([[self.rates[i,j]]*np.floor(self.fractions[i,j]*self.numLoci) for j in range(len(self.rates[i,]))])
        if (len(probRec)<self.numLoci):
            probRec=np.hstack((probRec,[self.rates[i,4]]*(self.numLoci-len(probRec))))
        elif (len(probRec)>self.numLoci):
            probRec=probRec[0:self.numLoci]
        probRec=probRec/sum(probRec)
        return probRec


class haploidGenome:
    '''
    haploid genome with a given number of loci(numLoci), mutation rate (mutRate), 
    and recombination probability at each locus (recProb). 
    Optional run length encoding for compression with useRLE.
    Currently only supports a single chromosome.
    '''
    def __init__(self, numLoci, mutRate, recProb, useRLE=False):
        self.numLoci=numLoci
        self.mutRate=mutRate
        self.recProb=recProb
        self.useRLE=useRLE
        self.loci=[]

    def getGenome(self):
        if (self.useRLE==True):
            return decode(self.loci)
        return self.loci

    def setGenome(self, genome):
        if (self.useRLE==True):
            self.loci=encode(genome)
        else:
            self.loci=genome.astype(bool)

    def mutate(self):
        mut=np.random.binomial(1,self.mutRate,self.numLoci)
        self.setGenome(abs(self.getGenome()-mut))

    def qtlGenotypeVal(self,qtlLoci,effectSize):
        return sum(self.getGenome()[qtlLoci]*effectSize)

class founderGenome(haploidGenome):
    '''
    a haploid genome of a parental stain.
    contains method to create a founder genome of a single genotype (0 or 1)
    of a specified length
    '''
    def __init__(self,haploidGenome):
        self.numLoci=haploidGenome.numLoci
        self.mutRate=haploidGenome.mutRate
        self.recProb=haploidGenome.recProb
        self.useRLE=haploidGenome.useRLE

    def createFounder(self, genotype):
        if (genotype==1):
            self.setGenome(np.ones(self.numLoci,dtype=bool))
        else:
            self.setGenome(np.zeros(self.numLoci,dtype=bool))
      

class daughterGenome(haploidGenome):
    '''
    a hapliod genome that inherits its attributes from parents. 
    contains method for recombining two parental genomes to generate a new 
    daughter genome
    '''
    def __init__(self,haploidGenome):
        self.numLoci=haploidGenome.numLoci
        self.mutRate=haploidGenome.mutRate
        self.recProb=haploidGenome.recProb
        self.useRLE=haploidGenome.useRLE
    
    def recombineOnce(self,genomeA,genomeB,mutate=True):
        if (mutate==True):
            genomeA.mutate()
            genomeB.mutate()
        breakpoint=np.random.choice(range(0,genomeA.numLoci),p=genomeA.recProb) 
        if (rnd.randint(0,1)==0):
            gamete=np.hstack((genomeA.getGenome()[0:breakpoint],
                           genomeB.getGenome()[breakpoint:]))
        else:
            gamete=np.hstack((genomeB.getGenome()[0:breakpoint],
                           genomeA.getGenome()[breakpoint:]))
        self.setGenome(gamete)

class Population:
    '''
    a collection of genomes
    '''
    def __init__(self, N, currentGen=0):
        self.N=N
        self.currentGen=currentGen
        self.genomes=[]

    def createF2(self,founder1,founder2,mutRate,recProb):
        F2=[]
        for g in range(0,self.N*2):
            gamete=daughterGenome(founder1)
            gamete.recombineOnce(founder1,founder2)
            F2.append(gamete)
        self.currentGen=2
        self.genomes=F2

    def createRIpop(self,toGeneration):
        if self.currentGen>=toGeneration:
            print("Population has already undergone %d generations" % 
                    self.currentGen)
        for gen in range(self.currentGen+1,toGeneration+1):
            tempPop=[]
            gameteIndex1=0
            gameteIndex2=1
            for indiv in range(self.N):
                genomeA=self.genomes[gameteIndex1]
                genomeB=self.genomes[gameteIndex2]
                gamete1=daughterGenome(genomeA)
                gamete2=daughterGenome(genomeA)
                gamete1.recombineOnce(genomeA,genomeB)                
                gamete2.recombineOnce(genomeA,genomeB)
                tempPop.append(gamete1)
                tempPop.append(gamete2)
                gameteIndex1+=2
                gameteIndex2+=2
            self.genomes=tempPop
            rnd.shuffle(self.genomes)
        self.currentGen=toGeneration

    def pheSelect(self, qtlLoci, effectSize, selectionDepth):
        '''
        selects selectionDepth% individuals from population according to
        genotype value at qtl loci ("selected population")
        '''
        selected=selectedPopulation(self)
        selected.N=int(self.N*selectionDepth/100)
        genotypeVals=np.array([g.qtlGenotypeVal(qtlLoci,effectSize) 
        for g in self.genomes],dtype=int).reshape(self.N,2).transpose()
        diploidVals=genotypeVals[0]+genotypeVals[1]
        varEnv=np.var(diploidVals)
        pheVals=diploidVals+np.random.normal(0,varEnv,self.N)
        pheOrder=pheVals.argsort()[::-1][:selected.N]
        indexRow=list(pheOrder*2)+list(pheOrder*2+1)
        selected.genomes=[self.genomes[i] for i in indexRow]
        return selected 
    
    def randomSelect(self,selectionDepth):
        '''
        selects selectionDepth% individuals from population randomly ("control
        population")
        '''
        selected=selectedPopulation(self)
        selected.N=int(self.N*selectionDepth/100)
        index=np.array(rnd.sample(range(self.N),selected.N))
        indexList=list(index*2)+list(index*2+1)
        selected.genomes=[self.genomes[i] for i in indexList]
        return selected


    def seqSample(self,seqDepth):
        '''
        randomly picks genomes at each locus independantly. this is a better 
        approximation of sequencing a DNA sample. always use sampling with 
        replacement as it is a better approximation of sequencing from pool of 
        DNA and avoids strange artefacts in different populations sizes
        '''
        numChrs=len(self.genomes)
        numLoci=self.genomes[0].numLoci
        freq=[]
        for j in range(numLoci):
            toSeq=[self.genomes[i] for i in np.random.choice(numChrs,seqDepth)]
            freq.append(np.array([g.getGenome()[j] for g in toSeq],dtype=int).mean())
        return np.array(freq)
        

class selectedPopulation(Population):
    '''
    a subset of genomes from a preexisting population
    '''
    def __init__(self,Population):
        self.N=0
        self.currentGen=Population.currentGen
        self.genomes=[]

def writeToTxtFile(contFreq,selectFreq,baseDir,attributes,values):
    contfname=baseDir+'cont_'+attributes.format(*values)
    selectfname=baseDir+'select_'+attributes.format(*values)
    with open(contfname,'at') as f:                   
        f.write(' '.join(['%1.6g' %i for i in contFreq])+'\n')
    with open(selectfname,'at') as f:                   
        f.write(' '.join(['%1.6g' %i for i in selectFreq])+'\n')
        
def writePopToFile(population,baseDir,numSim):
    fullDir=baseDir+'pop_N'+str(population.N)+'_G'+str(population.currentGen)+'/'
    try:
        os.stat(fullDir)
    except:
        os.mkdir(fullDir) 
    popfile=fullDir+'sim'+str(numSim)+'.npz'
    np.savez(popfile,*[g.loci for g in population.genomes])

def readPopFromFile(filename,haploidGenome):
    npzfiles=np.load(filename)
    fields=filename.split("/")
    attributes=fields[-2].split("_")
    N=int(attributes[1][1:])
    G=int(attributes[2][1:])
    population=Population(N,currentGen=G)
    for f in npzfiles.files:
        genome=daughterGenome(haploidGenome)
        genome.setGenome(decode(npzfiles[f]))
        population.genomes.append(genome)
    return population        
        