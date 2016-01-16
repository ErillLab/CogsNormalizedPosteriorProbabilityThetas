# -*- coding: utf-8 -*-
#imports the NCBI_utils lib containing the functions for GI accession and BLAST
import MG_synth_lib as MGlib

#random number generation
import random
#mean and std
from numpy import mean
from numpy import std
#normal dist
from scipy.stats import norm
#system calls
import sys

#set default parameters
motif_filename="CsoR.txt"   #input file
out_filename="_o"           #o prefix for output
verbose=0                   #verbose mode
alpha=1.0/300.0             #mixing ratio for regulated model
rproms=3.0                  #number of regulated promoters [prior]
tproms=1811.0               #total number of promoters in genome [prior]
experiment=3                #the experiment number

out_filename=motif_filename.split(".")[0] + out_filename + str(experiment)
    
#verbose
if verbose: print "Using: ", motif_filename, " as input"
if verbose: print "Writing to (suffix): ", "[void]" if out_filename==""\
else out_filename

#open files for ouput
try:
    out_file = open(out_filename + ".csv","w")
except (IOError, OSError) as file_open_exception:
    print "*** Something went wrong while opening the output file"
    print "*** Error: ", file_open_exception.errno, " - ",\
                         file_open_exception.strerror        
    sys.exit()      
    
#compute priors
PR=rproms/tproms               #prior probability of regulation
PB=1.0-PR                      #prior probability of non-regulation
PPR=PB/PR                      #prior probability ratio
    
#read motif and assing 0.25 pseudocounts to PSWM
#also assign background uniform distribution for the PSSM (default)
mot = MGlib.read_motif(motif_filename)
mot.pseudocounts=1
mot.background=None

#save the pssm for the motif and the reverse complement
#(so that they are not recalculated everytime we invoke motif.pssm)
pssm = mot.pssm
rpssm = pssm.reverse_complement()

#-------------------------------------------------------------------------
#Experiment 3:
#Using 100 sequences, with 3, 6, 9, 12, 15, 18 sites inserted
#get unormalizazed probability
#repeat 1000 times
random.seed(None)
verbose=1

#write csv header
out_file.write('Ins_sites,Post\n')

#loop experiments
for cnt in range(0,1000):
    print "Experiment: ", cnt

    #create background sequence set: 100 seqs 283 bp long
    set1 = MGlib.random_DNA(300,{'A': 0.3,'C': 0.2,'G': 0.2,'T': 0.3},100)
    set2 = set1[:]
    #compute softmax scores for background sequences in dataset
    gscr = MGlib.esfmax_score_seqs(set1,pssm,rpssm)
    #compute softmax scores for motif sequences
    mscr = MGlib.esfmax_score_seqs(mot.instances,pssm,rpssm)
    
    #get normal distributions for background and motif
    n_g=norm(mean(gscr), std(gscr))
    n_m=norm(mean(mscr), std(mscr))
   
    #create motif instances
    pmot1 = MGlib.sample_motif(mot,100)

    #insert sites in sequences
    e=0
    while (e<len(set1)):
        #edit sequence to include random site 
        set1[e] = random.choice(pmot1) + set1[e][17:]
        e = e+1
        
    #compute softmax scores for sequences in both datasets
    scrs1=MGlib.esfmax_score_seqs(set1,pssm,rpssm)
    scrs2=MGlib.esfmax_score_seqs(set2,pssm,rpssm)

    #in each experiment, create dataset with 100 seqs, sz with sites inserted
    for cnt2 in range(0,6):
        #determine dataset size
        if (cnt2==0): 
            sz=3
        elif (cnt2==1):
		sz=6
        elif (cnt2==2):
            sz=9
        elif (cnt2==3):
            sz=12
        elif (cnt2==4):
            sz=15
        else:
            sz=18
            
        #create mixed dataset score
        curr_scores=scrs1[0:sz]+scrs2[sz:len(scrs1)]

        #get log-likelihoods for sequences in dataset        
        llrs1=MGlib.ll_ratios(curr_scores,n_g,n_m,alpha)
        
        #get overall normalized posterior for the sequences in dataset
        posts1=MGlib.PostP(llrs1,PPR)

        if verbose: print sz, " - ", posts1

        #write results to file
        out_file.write(str(sz))
        out_file.write(',')
        out_file.write(str(posts1))
        out_file.write('\n')
            
        
out_file.close()