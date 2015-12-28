# BSA_simulation
Simulation of bulk segregant analysis - control vs selected large pools of recombinant inbred individuals

Acknowledgements: I initially developed simulations using the Hypred software in R by Frank Technow (http://finzi.psych.upenn.edu/library/hypred/html/hypred-package.html), but when it was not supported with updated versions of R I re-wrote the simulations in python to be faster. The code is all my own, but some of the conceptual thinking about how to do the simulation is drawn from Hypred.

Simulations are done in two stages: 
1) recombinant populations are created with runCreatePop.py using the functions defined in recombSim2.py. You can edit the script to specify the size of the population (Ns), number of simulations (numSims) and values for the number of generations (Gs) at which you want to store the populations generated. These generations are stored in compressed numpy format (.npz) in separate folders by population size and generations number. 

With populations larger than 10,000 it is best to run simulations in parallel on a cluster as they can take many hours to run (100,000 runs in almost 20 hours per simulation)

2) To test specific questions about experimental parameters, the recombinant populations (.npz files) are used to simulate phenotypic selection and sequencing to generate allele frequencies at each locus in each experiment, using the runSelect_XXXX.py scripts. I have created separate versions to test the effect of population size (runSelect_popn.py), number of generations (runSelect_generation.py), % of population selected (runSelect_selectDepth.py) and sequencing depth (runSelect_seqDepth.py). These scripts also require the funcitons in recombSim2.py.
The output is a set of pairs o text files (one for the control population and one for the selected population) for each combination of parameters (that are recorded in the filename). Each row in the file contains allele frequencies at each locus in a single simulation (50 sims means 50 rows). 
