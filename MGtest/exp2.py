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
experiment=2                #the experiment number

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
#Experiment 2:
#Using 100 sequences, with 12 sites inserted
#get normalizazed probability
#do this for 2, 3, 4, 5, 6 and 8 stdevs below mean
#repeat 1000 times
random.seed(None)
verbose=1

#write csv header
out_file.write('Theta,Ins_sites,Post,Pass_filt_seqs,Def_post\n')

#loop experiments
for cnt in range(0,1000):
    print "Experiment: ", cnt
    #in each experiment, create dataset with 100 seqs, sz with sites inserted
    for cnt2 in range(0,1):
        #create background sequence set: 100 seqs 283 bp long
        set1 = MGlib.random_DNA(283,{'A': 0.3,'C': 0.2,'G': 0.2,'T': 0.3},100)
        #compute softmax scores for background sequences in dataset
        gscr = MGlib.esfmax_score_seqs(set1,pssm,rpssm)
        #compute softmax scores for motif sequences
        mscr = MGlib.esfmax_score_seqs(mot.instances,pssm,rpssm)
        
        #get normal distributions for background and motif
        n_g=norm(mean(gscr), std(gscr))
        n_m=norm(mean(mscr), std(mscr))
       
        #create motif instances
        pmot1 = MGlib.sample_motif(mot,250)
        
        #determine dataset size
        if (cnt2==0): 
            sz=12
            
        #insert sites in sequences
        e=0
        while (e<len(set1)):
            #insert random site in first 10 sequences 
            if (e<sz+1): 
                set1[e] = random.choice(pmot1) + set1[e]
            #otherwise insert random sequence from own sequence start
            else :
                set1[e] = set1[e][:17] + set1[e]
            e = e+1

        #compute softmax scores for sequences in dataset
        scrs1=MGlib.esfmax_score_seqs(set1,pssm,rpssm)

        #perform calculations for each theta value
        for cnt3 in range(0,6):
            if (cnt3==0):
                theta=-8
            elif (cnt3==1):
                theta=-6
            elif (cnt3==2):
                theta=-5
            elif (cnt3==3):
                theta=-4
            elif (cnt3==4):
                theta=-3
            else:
                theta=-2
                
            #NORMALIZATION
            #compute effective cutoff (th) value given theta
            th = MGlib.ThetaTh(theta,n_m)
            if verbose: print "Effective cut-off: ", th                
            
            #compute revised priors (assuming 300 bp average length)
            aPR, aPB = MGlib.NormPriors(th, n_g, n_m, alpha, rproms,\
                                        tproms, promlen=300.0)
            #get revised prior ratio
            aPPR = aPB/aPR


            #FILTER sequences not matching theta
            #get list of sequences with min score > th
            Nscrs1 = [x for x in scrs1 if max(x)>=th]

            if verbose: print "Length of post-theta bckg seqs: ", len(Nscrs1)

            #get log-likelihoods for sequences in dataset        
            Nllrs1=MGlib.ll_ratios(Nscrs1,n_g,n_m,alpha)
            
            #get normalization factors
            Nnormfs1=MGlib.lNormFactor(Nscrs1, th, n_g, n_m, alpha)

            #get overall normalized posterior for the sequences in dataset
            Nposts1=MGlib.NormPostP(Nllrs1,aPPR,Nnormfs1,0)

            if verbose: print theta, " - ", sz, " - ", Nposts1, " - ", \
                              len(Nscrs1), 1/(1+aPPR)
            #write results to file
            out_file.write(str(theta))
            out_file.write(',')
            out_file.write(str(sz))
            out_file.write(',')
            out_file.write(str(Nposts1))
            out_file.write(',')
            out_file.write(str(len(Nscrs1)))
            out_file.write(',')
            out_file.write(str(1/(1+aPPR)))
            out_file.write('\n')
            
        
out_file.close()