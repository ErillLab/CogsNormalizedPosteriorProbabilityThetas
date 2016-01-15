# -*- coding: utf-8 -*-
"""
This class provides a simple wrapper for Biopython's PSSM scoring functions.

Created on Wed Feb 11 11:05:08 2015

@author: Talmo
"""

import os
import random
from Bio import Seq, SeqRecord, Alphabet, motifs
import numpy as np
import scipy
import scipy.stats


class PSSMScorer:
    def __init__(self, binding_sites, name="", pseudocounts=1):
        """ Creates a PSSM scorer object. Accepts binding sites in the form of
        a path to a text file containing one site per line, or a list of 
        Biopython Bio.Seq objects. """
        
        self.name = name
        self.alphabet = Alphabet.IUPAC.unambiguous_dna
        self.path = None
        
        if type(binding_sites) == str:
            self.seqs = [Seq.Seq(site.strip(), self.alphabet) for site in open(binding_sites).readlines()]
            self.name = os.path.splitext(os.path.basename(binding_sites))[0]
            self.path = binding_sites
        elif type(binding_sites) == list:
            self.seqs = binding_sites
            
        self.n = len(self.seqs)
        
        # Default name
        if len(self.name) == 0:
            self.name = "pssm_%dbp_%dseqs" % (self.m, len(self.seqs))
            
        # Construct motif
        self.motif = motifs.create(self.seqs)
        self.motif.pseudocounts = pseudocounts
        
        # Construct PSSM and reverse PSSM
        self.pssm = self.motif.pssm
        self.pssm_r = self.pssm.reverse_complement()
        self.m = self.pssm.length
        self.w = self.pssm.length
        self.length = self.pssm.length
        
        # Fast score primitives
        self.dict_pssm = dict(self.pssm)
        self.dict_pssm_r = dict(self.pssm_r)
        
        # Bayesian estimator
        self.estimator_initialized = False
        
    def __repr__(self):
        return "%s [%d bp | %d seqs]" % (self.name, self.length, self.n)
    def __str__(self):
        return "%s [%d bp | %d seqs]" % (self.name, self.length, self.n)
    def __iter__(self):
        return iter(self.seqs)
        
    def convert_seq(self, seq):
        """ Converts sequence strings to Biopython Seq objects with appropriate
        alphabet. """
        
        if type(seq) == str:
            return Seq.Seq(seq.strip(), self.alphabet)
        elif type(seq) == SeqRecord.SeqRecord:
            return self.convert_seq(seq.seq)
        elif type(seq) == Seq.Seq and seq.alphabet != self.alphabet:
            seq.alphabet = self.alphabet
            return seq
        else:
            return seq
        
    def score_bio(self, seq):
        """ Scores a sequence using Biopython and returns the best score between the forward
        and reverse strand. """
        if len(seq) != self.m:
            raise Exception("Sequence must be of same length as PSSM.")
            
        seq = self.convert_seq(seq)
        return max(self.pssm.calculate(seq), self.pssm_r.calculate(seq))
        
    def score_all(self, seq):
        """ Scores all sites in a sequence and returns an array of scores. """
        all_scores = self.search(seq, -np.inf)
        scores = all_scores[all_scores[:, 0] >= 0, 1]
        scores_r = all_scores[all_scores[:, 0] < 0, 1]
        
        return scores, scores_r
        
    def search(self, seq, threshold=0.0):
        """ Search for the sites in the sequence with a score above a threshold.
        Searches on both strands."""
        return np.array(list(self.pssm.search(self.convert_seq(seq), both=True, threshold=threshold)))

    def score(self, seq, soft_max=False):
        """ Sliding window scorer using fast primitives. """
        n = len(seq)
        scores = np.zeros(n - self.m + 1)
        scores_r = np.zeros(n - self.m + 1)
        for i in xrange(n - self.m + 1):
            for pos in xrange(self.m):
                scores[i] += self.dict_pssm[seq[i+pos]][pos]
                scores_r[i] += self.dict_pssm_r[seq[i+pos]][pos]
        if soft_max:
            return sm(scores, scores_r)
        else:
            return scores, scores_r

    def score_self(self, soft_max=False):
        """ Scores the sequences used to build the motif. """
        scores = [self.score(seq, soft_max) for seq in self.seqs]
        
        if soft_max:
            scores = np.hstack(scores) # convert to vector
        return scores

    def initialize_estimator(self, bg_mu=None, bg_sigma=None, num_random=100000):
        """ Initializes the parameters for the Bayesian estimator. """
        # Parameters
        self.pf = 1 / 100.0 # foreground probability
        self.pb = 1 - self.pf # background probability
        self.alpha = 1.0 / 300 # frequency of binding site?
        #num_random = 10000 # number of random sequences to score to estimate background
        
        # Score motif sequences to estimate foreground
        pssm_scores = self.score_self(True)
        mu_y, sigma_y = np.mean(pssm_scores), np.std(pssm_scores)
        
        # Background
        mu_x = bg_mu
        sigma_x = bg_sigma
        if bg_mu is None or bg_sigma is None:
            # Generate random background scores
            # TODO: Read this off the PSSM
            background_scores = np.array([self.score(random_seq(self.length), True) for i in xrange(int(num_random))])
            if bg_mu is None:
                mu_x = np.mean(background_scores)
            if bg_sigma is None:
                sigma_x = np.std(background_scores)
        
        # Distributions
        self.pdf_y = scipy.stats.distributions.norm(mu_y, sigma_y).pdf
        self.pdf_x = scipy.stats.distributions.norm(mu_x, sigma_x).pdf
    
        # Calculations
        self.L_b = lambda scores: self.pdf_x(scores)
        self.L_f = lambda scores: self.alpha * self.pdf_y(scores) + (1 - self.alpha) * self.pdf_x(scores)
        self.LL_f = lambda scores: np.exp(np.sum(np.log(self.L_f(scores))))
        self.LL_b = lambda scores: np.exp(np.sum(np.log(self.L_b(scores))))
        self.LL_ratio = lambda scores: np.exp(np.sum(np.log(self.L_b(scores)) - np.log(self.L_f(scores))))
        
        # Update initialized flag
        self.estimator_initialized = True
        
    def post_prob(self, sm_scores):
        """ Computes the posterior probability that the scores contain a
        binding site. """
        if not self.estimator_initialized:
            self.initialize_estimator()
            
        return 1 / (1 + self.LL_ratio(sm_scores) * self.pb / self.pf)
        
#%% Static methods
def random_seq(l):
    """ Generates a random sequence sampling from the uniform distribution. """
    return "".join([random.choice("ACGT") for i in xrange(l)])
        
def sm(scores, scores_r):
    """ Computes the soft max of the scores in two strands. """
    return np.log(np.exp(scores) + np.exp(scores_r))
    
