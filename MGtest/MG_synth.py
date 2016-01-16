# -*- coding: utf-8 -*-
"""
Created on Wed Dec 2 2015

Main for MG_Synth. This main calls the necessary functions to create a
synthetic dataset for analysis of the metagenome regulatory network Bayesian 
inference module

@author: ivanerill
"""
#imports getopt to handle the cmd line reading
import getopt
#imports the NCBI_utils lib containing the functions for GI accession and BLAST
import MG_synth_lib as MGlib

#system calls
import sys
#random number generation
import random
#mean and std
from numpy import mean
from numpy import std
#normal dist
from scipy.stats import norm

def main():
    """Gets a motif from file and reads it. It then generates synthetic data
       to represent a set of promoters (100) mapping to a particular eggNOG/COG,
       inserts into these sequences pseudosites (generated from the 
       distribution implicit in the PSSM). It then calls the PSSM evaluation
       function to score the sites using the softmax function and then the
       different functions to compute the likelihoods and the posterior
       probabilities.
       
       Usage:
       
       MG_synth -M <Motif file> -O <out file prefix> -E <experiment> \
                -A <alpha mix ratio> -P <Regulation prior> -T <theta> \
                -V <verbose mode>
       
       
       Note: motifs are assumed to be in FASTA or 1-per-line text format
    """
     
    #set default parameters
    motif_filename="CsoR.txt"   #input file
    out_filename="_o"           #o prefix for output
    verbose=0                   #verbose mode
    alpha=1.0/300.0             #mixing ratio for regulated model
    rproms=3.0                  #number of regulated promoters [prior]
    tproms=1811.0               #total number of promoters in genome [prior]
    experiment=2                #the experiment number
    
    #get cmd parameters
    try:
        opts, args=getopt.getopt(sys.argv[1:],"I:O:V")
    except getopt.GetoptError:
        print 'MG_synth -M <Motif file> -O <out file prefix> -E <experiment> \
                -A <alpha mix ratio> -P <Regulation prior> -T <theta> \
                -V <verbose mode>'
        sys.exit(2) 
        
    #assign parameters
    for opt, arg in opts:
        if opt == '-M':
            motif_filename=arg
        elif opt == '-O':
            out_filename=arg
        elif opt == '-E':
            experiment=int(arg)
        elif opt == '-A':
            alpha=float(arg)
        elif opt == '-P':
            PR=float(arg)
        elif opt == '-T':
            theta=float(arg)
        elif opt == '-V':
            verbose=arg
        elif opt == '-askme':
            motif_filename = raw_input('Enter the motif file name\n')
            out_filename = raw_input('Enter the output file name prefix\n')
            experiment = raw_input('Enter the experiment number\n')
            alpha = raw_input('Enter the alpha mixing ratio\n')
            PR = raw_input('Enter the prior probability for regulation\n')
            theta = raw_input('Enter the theta sensitivity threshold\n')
            verbose = raw_input('Enter verbose mode (1/0)\n')
            
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
    
    #open file for error recording
    try:
        err_file = open(out_filename+".err","w")
    except (IOError, OSError) as file_open_exception:
        print "Something went wrong while opening the error file"
        sys.exit()

    #compute priors
    PR=rproms/tproms               #prior probability of regulation
    PB=1.0-PR                      #prior probability of non-regulation
    PPR=PB/PR                      #prior probability ratio
        
    #read motif and assing 0.25 pseudocounts to PSWM
    #also assign background uniform distribution for the PSSM (default)
    mot = MGlib.read_motif(motif_filename)
    mot.pseudocounts=0.25
    mot.background=None
    
    #save the pssm for the motif and the reverse complement
    #(so that they are not recalculated everytime we invoke motif.pssm)
    pssm = mot.pssm
    rpssm = pssm.reverse_complement()

    if (experiment==0):
        #-------------------------------------------------------------------------
        #Experiment 0:
        #10000 sequences, with 100% on average having random pseudo-sites inserted
        #all sites are inserted at the first position of the sequence
        #get individual sequence posteriors and write them together with the
        #score of the site inserted
        random.seed(None)
        
        #create background sequence set: 100 seqs 283 bp long
        set1 = MGlib.random_DNA(283,{'A': 0.3,'C': 0.2,'G': 0.2,'T': 0.3},10000)
        #compute softmax scores for background sequences in dataset
        gscr = MGlib.sfmax_score_seqs(set1,pssm,rpssm)
     
        #get normal distributions for background and motif
        n_g=norm(mean(gscr), std(gscr))
        n_m=norm(pssm.mean(), pssm.std())
       
        #create motif instances
        pmot1 = MGlib.sample_motif(mot,1000)

        #insert sites in sequences
        e=0
        while (e<len(set1)):
            r = random.random()
            #determine if site is to be inserted and insert random site
            if (r<1): 
                set1[e] = random.choice(pmot1) + set1[e]
            #otherwise insert random sequence from own sequence start
            else :
                set1[e] = set1[e][:17] + set1[e]
            e = e+1

        #compute softmax scores for sequences in dataset
        scrs=MGlib.sfmax_score_seqs(set1,pssm,rpssm)
        #get log-likelihoods for sequences in dataset        
        llrs=MGlib.ll_ratios(scrs,n_g,n_m,alpha)
        
        #get per-sequence posterior for the sequences in dataset
        fposts=MGlib.PostP(llrs,PPR,1)

        #write results to file
        e=0
        while (e<len(set1)):
            out_file.write(str(scrs[e][0]))
            out_file.write(',')
            out_file.write(str(fposts[e]))
            out_file.write('\n')
            e=e+1

        out_file.close()
        return 0
    
    elif (experiment==1):
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
        
        for cnt in range(0,100):
            #create background sequence set: 100 seqs 283 bp long
            set1 = MGlib.random_DNA(283,{'A': 0.3,'C': 0.2,'G': 0.2,'T': 0.3},100)
            set2 = set1[:]
            #compute softmax scores for background sequences in dataset
            gscr = MGlib.sfmax_score_seqs(set1,pssm,rpssm)
         
            #get normal distributions for background and motif
            n_g=norm(mean(gscr), std(gscr))
            n_m=norm(pssm.mean(), pssm.std())
           
            #create motif instances
            pmot1 = MGlib.sample_motif(mot,100)
            pmot2 = MGlib.sample_motif(mot,1)
            
            print cnt
            #insert sites in sequences
            e=0
            while (e<len(set1)):
                r = random.random()
                #insert random site in first 10 sequences 
                if (e<11): 
                    set1[e] = random.choice(pmot1) + set1[e]
                    set2[e] = random.choice(pmot2) + set2[e]
                #otherwise insert random sequence from own sequence start
                else :
                    set1[e] = set1[e][:17] + set1[e]
                    set2[e] = set2[e][:17] + set2[e]
                e = e+1
    
            #compute softmax scores for sequences in dataset
            scrs1=MGlib.sfmax_score_seqs(set1,pssm,rpssm)
            scrs2=MGlib.sfmax_score_seqs(set2,pssm,rpssm)
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
        return 0
        
   elif (experiment==2):
        #-------------------------------------------------------------------------
        #Experiment 2:
        #Using 100 sequences, with 2, 2, 5, and 10% of sites inserted
        #get unormalizazed and normalized probability
        #do this for 1, 2, 4 and 6 stdevs below mean
        #repeat 100 times
        random.seed(None)
        verbose=1
        print "Experiment #2"
        for cnt in range(0,1000):
            for cnt2 in range(0,4):
                #create background sequence set: 100 seqs 283 bp long
                set1 = MGlib.random_DNA(283,{'A': 0.3,'C': 0.2,'G': 0.2,'T': 0.3},100)
                #compute softmax scores for background sequences in dataset
                gscr = MGlib.sfmax_score_seqs(set1,pssm,rpssm)
             
                #get normal distributions for background and motif
                n_g=norm(mean(gscr), std(gscr))
                n_m=norm(pssm.mean(), pssm.std())
               
                #create motif instances
                pmot1 = MGlib.sample_motif(mot,250)
    
                print cnt
                
                #determine dataset size
                if (cnt2==0): 
                    sz=5
                elif (cnt2==1):
                    sz=10
                elif (cnt2==2):
                    sz=15
                elif (cnt2==3):
                    sz=20
                    
                #insert sites in sequences
                e=0
                while (e<len(set1)):
                    r = random.random()
                    #insert random site in first 10 sequences 
                    if (e<sz+1): 
                        set1[e] = random.choice(pmot1) + set1[e]
                    #otherwise insert random sequence from own sequence start
                    else :
                        set1[e] = set1[e][:17] + set1[e]
                    e = e+1
        
                #compute softmax scores for sequences in dataset
                scrs1=MGlib.sfmax_score_seqs(set1,pssm,rpssm)

                for cnt3 in range(0,4):
                    if (cnt3==0):
                        theta=-6
                    elif (cnt3==1):
                        theta=-4
                    elif (cnt3==2):
                        theta=-2
                    else:
                        theta=-1
                        
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
    
                    if verbose: print "Length of post-theta bckg seqs: ",\
                                len(Nscrs1)
    
                    #get log-likelihoods for sequences in dataset        
                    Nllrs1=MGlib.ll_ratios(Nscrs1,n_g,n_m,alpha)
                    
                    #get normalization factors
                    Nnormfs1=MGlib.lNormFactor(Nscrs1, th, n_g, n_m, alpha)
    
                    #overall normalized posterior for the sequences in dataset
                    Nposts1=MGlib.NormPostP(Nllrs1,aPPR,Nnormfs1,0)

                    if verbose: print theta, " - ", sz, " - ", Nposts1, " - ",\
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
        return 0
        
    
if __name__ == '__main__':
    main()