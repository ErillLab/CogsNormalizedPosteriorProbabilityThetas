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
experiment=1                #the experiment number

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
#Experiment 1:
#2x100 sequences, with the first 10 having random pseudo-sites inserted
#all sites are inserted at the first position of the sequence
#get individual sequence posteriors and write them together with the
#score of the site inserted
#first sequence set includes randomly distributed sites
#second sequence set includes only one site
#experiment is repeated 100 times
random.seed(None)

#write csv header
out_file.write('Post_motif,Post_single_site\n')

#loop experiments
for cnt in range(0,100):
    print "Experiment: ", cnt
    #create background sequence set: 100 seqs 283 bp long
    set1 = MGlib.random_DNA(283,{'A': 0.3,'C': 0.2,'G': 0.2,'T': 0.3},100)
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
    pmot2 = MGlib.sample_motif(mot,1)
    
    #insert sites in sequences
    e=0
    while (e<len(set1)):
        #insert random site in first 10 sequences 
        if (e<11): 
            set1[e] = random.choice(pmot1) + set1[e]
            set2[e] = pmot2[0] + set2[e]
        #otherwise insert random sequence from own sequence start
        else :
            set1[e] = set1[e][:17] + set1[e]
            set2[e] = set2[e][:17] + set2[e]
        e = e+1

    #compute softmax scores for sequences in dataset
    scrs1=MGlib.esfmax_score_seqs(set1,pssm,rpssm)
    scrs2=MGlib.esfmax_score_seqs(set2,pssm,rpssm)
    if verbose: print "varied"
    if verbose: 
        for s in scrs1: print s[0]
    if verbose: print "loners"
    if verbose: 
        for s in scrs2: print s[0]
    #get log-likelihoods for sequences in dataset        
    llrs1=MGlib.ll_ratios(scrs1,n_g,n_m,alpha)
    llrs2=MGlib.ll_ratios(scrs2,n_g,n_m,alpha)
    
    #get overall posterior for the sequences in dataset
    fposts1=MGlib.PostP(llrs1,PPR,0)
    fposts2=MGlib.PostP(llrs2,PPR,0)

    #write results to file
    out_file.write(str(fposts1))
    out_file.write(',')
    out_file.write(str(fposts2))
    out_file.write('\n')
    
out_file.close()

