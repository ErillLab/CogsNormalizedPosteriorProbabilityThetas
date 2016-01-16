# -*- coding: utf-8 -*-
"""
Created on Wed Dec 2 2015

Library for MG_Synth. Defines functions for reading motifs, generating 
synthetic DNA sequences, computing likelihoods, etc

@author: ivanerill
"""

#The following code assumes that Biopython is installed
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC
#system calls
import sys
#random number generation
import random
#math functions
import math

#------------------------------------------------------------------------------
def read_motif(motif_filename, verb=0):
    """Reads a motif as a collection of sites from a file

    Reads a motif and uses the biopython.motifs class to store it. If the motif
    is in FASTA format, it uses the parser directly. Otherwise, it loads and 
    reads a concatenated text file and creates the motif.
    
    File type is determined by extension: 
        * FASTA for .fas, .fasta and .fa files
        * One-per-line text file otherwise
    
    Input:
        * The motif filename; required
        * Verbose mode (default=0)
    Returns:
        * the read motif
    """


    #create file handler for reading file
    try:
        motif_file = open(motif_filename,"r")
    except (IOError, OSError) as file_open_exception:
        print "*** The file name provided:", motif_filename, " does not exist"
        print "*** Error: ", file_open_exception.errno, " - ",\
                             file_open_exception.strerror       
        sys.exit()

    #Figure out file type based on extension, read sites and create motif
    extension = motif_filename.split('.')[-1]

    if extension not in ['fas', 'fasta', 'fa']:
        if verb: print 'Reading motif... raw text sequence mode assumed \
                        (one site per line, not FASTA parsing)'
        sites = []
        for line in motif_file:
            sites.append(Seq(line.rstrip('\n\r'),IUPAC.unambiguous_dna))
        mot = motifs.create(sites)
        if verb: print mot.degenerate_consensus
    else:
        if verb: print 'Reading motif... attempting to parse FASTA file'
        mot = motifs.read(motif_file,'sites')
        
    motif_file.close()
    return mot

#------------------------------------------------------------------------------
def sample_motif(mot, sampleN, verb=0):
    """
    Samples sampleN instances from a motif, using the positional-frequencies
    of the motif to generate proportional random instances of the motif
    
    Input:
        * the motif; required
        * the number of samples to be generated; required
        * Verbose mode (default=0)
    Returns:
        * the samples generated, as a list
    """

    random.seed(None)
    
    #the list to return
    samples = []
    #the column holder
    cols = []
    
    ind=0
    #for each position in the motif
    while ind<len(mot):
        #get the frequency list (dict) for the column and sort descending
        sorted_freqs=sorted(mot.pwm[:,ind].items(),key=lambda x: x[1],\
                            reverse=True)
        cnt = 0    
        column = []
        #for sample number to be drawn
        while cnt < sampleN:
            r  = random.random()
            cumsum = 0
            #determine base to draw according to uniform random and freq list
            for sf in sorted_freqs:
                cumsum = cumsum + sf[1]
                if (r<cumsum):
                    base = sf[0]
                    break
            #append the base to the column
            column.append(base)
                    
            cnt = cnt + 1
            
        
        #append column to column list
        #cols is a list of len(mot) elements, with each element 
        #being a list containing sampleN bases
        #ex: [['A','A','G'],['G','T','T']] could be for sampleN=3 and length=2
        #where a first position is dominated by A and the second by T
        if verb: print column
        cols.append(column)
        
        ind=ind+1
    
    #zip appended columns
    #this will unpack (*) the elements of cols: ['A','A','G'] and ['G','T','T']
    #and then zip them (which will merge the i-th element of each list and 
    #return it as a tuple): [('A', 'G'), ('A', 'T'), ('G', 'T')]
    cols=zip(*cols)
    
    #create samples
    #this iterates through the list of tuples return by zip, and merges
    #them into a string: ('A', 'G') becomes 'AG', the first motif instance
    #'AT' the second, and 'GT' the third
    for c in cols:
        samples.append(Seq(''.join(c),IUPAC.unambiguous_dna))
        
    #return samples
    #['AG','AT','GT']
    return samples
    
#------------------------------------------------------------------------------
def random_DNA(length, freqs, N, verb=0):
    """
    Generates N random DNA sequences using the provided %GC content
    
    Input:
        * the length of the sequence to be generated; required
        * the frequencies of each base as dictionary [A C G T]; required
          ex: {'A': 0.96, 'C': 0.0, 'G': 0.0, 'T': 0.04}
        * the number of sequences to be generated
        * Verbose mode (default=0)
    Returns:
        * the generated sequences as a list of biopython seq objects
    """

    random.seed(None)

    seqlist=[]
    sorted_freqs=sorted(freqs.items(),key=lambda x: x[1],reverse=True)
    
    #notify user if freqs do not add to 1
    if (sum(freqs.values())<>1):
        print "Frequencies do not add up to one!"
       
    ind=0
    #for each sequence to be generated
    while (ind<N):
        #for each base to be added
        cnt=0
        seq=''
        while (cnt<length):
            r = random.random()
            cumsum = 0
            for sf in sorted_freqs:
                cumsum = cumsum + sf[1]
                if (r<cumsum):
                    base = sf[0]
                    break
            #concatenate chosen base to growing DNA sequence
            seq=seq+base
            cnt=cnt + 1
            
        #make sequence object
        dnaseq=Seq(seq,IUPAC.unambiguous_dna)
        #append sequence object to list
        seqlist.append(dnaseq)
        ind = ind + 1
    
    return seqlist

#------------------------------------------------------------------------------
def sfmax_score_seqs(seq_list, pssm, rpssm=None, verb=0, mode=1):
    """
    Scores a list of sequences using provided pssm and reverse pssm, applies
    softmax function to return a single score per position
    
    Input:
        * the sequence list; required
        * the pssm to score with; required
        * the reverse pssm; computed if not provided
        * Verbose mode (default=0)

    Returns:
        * the softmax scores, as a list (one list of scores per sequence)
    """
    #list of list of scores to return
    scorelist=[]
    
    #compute rpssm if not provided
    if (rpssm==None): rpssm=pssm.reverse_complement()
    
    #compute scores for each sequence
    for s in seq_list:
        sc = pssm.calculate(s)
        rsc = rpssm.calculate(s)
        #handle sequences of length=motif length (0-dim array)
        if (sc.size==1):
            sc=[sc]
            rsc=[rsc]
        #apply softmax
        scores=map(lambda x: math.log(2.0**x[0]+2.0**x[1],2),zip(sc,rsc))
        #append list of scores to overall list (for each sequence)
        scorelist.append(scores)
    
    return scorelist

#------------------------------------------------------------------------------
def esfmax_score_seqs(seq_list, pssm, rpssm=None, verb=0, mode=1):
    """
    Scores a list of sequences using provided pssm and reverse pssm, applies
    softmax function to return a single score per position
    Uses the natural log (Boltzmann) derivation of the softmax
    
    Input:
        * the sequence list; required
        * the pssm to score with; required
        * the reverse pssm; computed if not provided
        * Verbose mode (default=0)
        
    Returns:
        * the softmax scores, as a list (one list of scores per sequence)
    """
    #list of list of scores to return
    scorelist=[]
    
    #compute rpssm if not provided
    if (rpssm==None): rpssm=pssm.reverse_complement()
    
    #compute scores for each sequence
    for s in seq_list:
        sc = pssm.calculate(s)
        rsc = rpssm.calculate(s)
        #handle sequences of length=motif length (0-dim array)
        if (sc.size==1):
            sc=[sc]
            rsc=[rsc]
        #apply softmax
        scores=map(lambda x: math.log(math.exp(x[0])+math.exp(x[1])),zip(sc,rsc))
        #append list of scores to overall list (for each sequence)
        scorelist.append(scores)
    
    return scorelist
#------------------------------------------------------------------------------
def ll_ratios(score_set, n_g, n_m, alpha, verb=0):
    """
    Computes the likelihood ratio for a given set of score lists (each a score 
    list corresponding to one sequence), using provided background and 
    regulated mean/stdevs and mixing ratio and assuming a normal distribution
    
    Input:
        * the score list; required
        * the background normal distribution; required
        * the regulated normal distribution; required
        * the mixing ratio; required
        * Verbose mode (default=0)
    Returns:
        * the log likelihood ratio, as a list (one per sequence) 
    """    
    
    #list of log-likelihood ratios to be returned, one per sequence/score list    
    llrs=[]
    
    #for each score list
    for score_list in score_set:
        #compute the sum of log likelihood ratios for score array
        sumlr=0.0
        for scr in score_list:
            lpd_b = n_g.logpdf(scr)
            lpd_f = math.log(alpha*n_m.pdf(scr) + (1-alpha)*n_g.pdf(scr))
            #add log-likelihood differences (i.e. mult. likelihood ratios)
            sumlr=sumlr + (lpd_b-lpd_f)
        llrs.append(sumlr)
    return (llrs)    

#------------------------------------------------------------------------------    
def PostP(LLR,PPR,mode=0,verb=0):
    """
    Returns the posterior probability for a set of sequences, given their
    log-likelihood ratios and the prior probability ratio
    
    Input:
        * the list of log-likelihood ratios; required
        * the prior probability ratio; required
        * mode: whether a list of posteriors or the aggregated is returned
        * Verbose mode (default=0)
    Returns:
        * the posterior probabilities, as a list (one per sequence)
    """
    
    if (mode):
        #list to be returned
        plist=[]    
        for loglr in LLR:
            plist.append(1.0/(1.0+math.exp(loglr)*PPR))
        return(plist)
    else:
        #return overall posterior
        cumllr=0
        for loglr in LLR:
            cumllr=cumllr+loglr
        return(1.0/(1.0+math.exp(cumllr)*PPR))

#------------------------------------------------------------------------------    
def ThetaTh(theta,n_m, verb=0):
    """
    Returns the sensitivity cut-off (th), given a theta value expressed as the
    number of standard deviations with respect to the mean for the distribution
    of known sites
    
    Input:
        * the minimum score, theta expressed as stdevs from mean motif score
        * the regulated normal distribution; required
        * Verbose mode (default=0)
    Returns:
        * the cutoff value
    """
    #compute effective cut-off score (th) as theta stdevs below mean
    mmean, mvar = n_m.stats(moments='mv') 
    th = mmean + theta * math.sqrt(mvar)
    
    return th

#------------------------------------------------------------------------------    
def NormPriors(th,n_g,n_m,alpha,rprom,tprom, promlen=300.0, verb=0):
    """
    Returns the normalized priors, given total and regulated number of promoters,
    the thresholding value th, the regulated and background models and an 
    assumed average length for promoter sequences
    
    The normalized priors are obtained as follows:
    - The number of non-regulated promoters is revised, given the expectation of 
      those promoters making it through the filtering step (i.e. we multiply
      the known prior times the probability of observing a sequence with score
      above th (after the filtering step) given the background model). This 
      gives us the expected number of total promoters in the new scenario
    - The number of regulated promoters is also revised, in accordance to the
      probability of observing a sequence with score above th given the mixture
      model. This gives us the expected number of regulated promoters in the new
      scenario.

    
    Input:
        * the minimum score, th
        * the background normal distribution; required
        * the regulated normal distribution; required
        * the mixing ratio; required
        * the number of regulated promoters; required
        * the total number of promoters; required
        * the average promoter sequence length; default=300
        * Verbose mode (default=0)
    Returns:
        * the normalized prior for regulation
    """
    
    #number of regulated promoters
    nrprom=tprom-rprom
    #get probability of observing a promoter with at least no score above
    #th, given the background model
    Ub = n_g.cdf(th) ** promlen
    Ur = (alpha*n_m.cdf(th)+(1.0-alpha)*n_g.cdf(th)) ** promlen
    
    if verb: print nrprom, " --> ", nrprom * (1-Ub)
    if verb: print rprom, " --> ", rprom * (1-Ur)
    
    #recompute the priors with the "expected" number of promoters after
    #discarding all promoters with no scores above th
    pr = (rprom * (1-Ur)) / ( (nrprom * (1-Ub)) + (rprom * (1-Ur)) )
    pnr = 1 - pr
    
    return [pr, pnr]
#------------------------------------------------------------------------------    
def lNormFactor(score_set, th, n_g, n_m, alpha, verb=0):
    """
    Returns, for a given set of scores, the log normalization ratios
    with sensitivity adjustment for some minimum score th taken under 
    consideration (one per sequence)
    These normalization ratios, which represent the log ratio between
    the probability of observing a sequence (score_list) with at least one
    score above the cutoff under the regulated and background models, will
    be used to adjust the likelihood ratios in the computation of the posterior
    
    Input:
        * the list of scores for each sequence
        * the minimum score, theta expressed as stdevs from mean motif score
        * the background normal distribution; required
        * the regulated normal distribution; required
        * the mixing ratio; required
        * Verbose mode (default=0)
    Returns:
        * the log normalization factors, as a list (one per sequence)
    """
    
    #the list of ratios to be returned
    lnormratios=[]

    #for each list of scores (i.e. sequence)
    for score_list in score_set:
        #compute Ub/Ur (prob. of observing score_list with no scores above th)
        #given background and regulated models
        Ub = n_g.cdf(th) ** len(score_list)
        Ur = (alpha*n_m.cdf(th)+(1.0-alpha)*n_g.cdf(th)) ** len(score_list)
        lnormratios.append(math.log((1.0-Ur)/(1.0-Ub)))
                
    return lnormratios
    
#------------------------------------------------------------------------------    
def NormPostP(LLR,PPR,lnormf,mode=0,verb=0):
    """
    Returns the posterior probability for a set of sequences, given their
    log-likelihood ratios, the prior probability ratio and their normalization
    factors
    
    Input:
        * the list of log-likelihood ratios; required
        * the prior probability ratio; required
        * the list of normalization factors; required
        * mode: whether a list of posteriors or the aggregated is returned
        * Verbose mode (default=0)
    Returns:
        * the posterior probabilities, as a list (one per sequence)
    """
    if (mode):
        #list to be returned
        plist=[]    
        for loglr, lnorm in zip(LLR, lnormf):
            plist.append(1.0/(1.0+math.exp(loglr+lnorm)*PPR))
        return(plist)
    else:
        #return overall posterior
        cumllr=0
        for loglr, lnorm in zip(LLR, lnormf):
            cumllr=cumllr+loglr+lnorm
        return(1.0/(1.0+math.exp(cumllr)*PPR))    


    
    
    
    
    
    