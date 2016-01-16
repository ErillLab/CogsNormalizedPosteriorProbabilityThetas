# MGtest
Library and test routines for metagenomics regulatory network Bayesian inference

##Library
The MG\_synth\_lib library contains the necessary routines to read a motif from 
file, sample instances from it, generate random DNA sequences according to a 
prespecified base composition, score sequences using the softmax function based 
on the motif PSSM, compute likelihood ratios for sequences and integrate them 
into posterior probabilities (given provided priors), and to normalize 
likelihood ratios, priors and posterior probabilities given a filtering step 
governed by the score threshold _theta_.

##Experiments
The experiments are as follows:

- Experiment 0: takes 100 randomly generated sequences with inserted sites and 
computes the individual sequence posterior probability for each, so that the 
relationship with inserted site score and posterior probability can be tracked
- Experiment 1: generates two sets of 100 sequences. It inserts randomly 
sampled instances of the motif in the first 10 sequences of one set, and a 
single sampled instance of the motif in the first 10 sequences of the other set.
It then computes the aggregated sequence set posteriors. This allows observing 
the dependence of the posterior on the distribution of sites.
- Experiment 2: generates 100 random sequences and inserts sites sampled from 
the motif in 12 of them. It computes the normalized probability for the sequence 
set under varying values of the score threshold _theta_ (2, 3, 4, 5, 6 and 8 
stdevs below the motif mean score). This is repeated 1000 times. This allows 
tracking the effect of _theta_ on the posterior probability. It also outputs 
the normalized prior for regulation and the number of filtered sequences.
- Experiment 3: generates 100 random sequences and inserts sites sampled from 
the motif in 3, 6, 9, 12, 15 and 18 of them. It computes the unnormalized 
probability for the sequence set. This is repeated 1000 times. This allows 
studying the dependence of the posterior on the relative number of regulated 
promoters within a sequence set.

