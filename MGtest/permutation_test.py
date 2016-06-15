# -*- coding: utf-8 -*-
#The following code assumes that Biopython is installed
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC

#imports the NCBI_utils lib containing the functions for GI accession and BLAST
import MG_synth_lib as MGlib

#random number generation
import random
#mean and std
from numpy import mean
from numpy import std
import numpy as np

#normal dist
from scipy.stats import norm
from scipy.stats import geom
#system calls
import sys
#import math # Using numpy instead for vectorized operations

###########################################################################################
def create_COG(mode, mot):
    # generate random sequence set1 (100 seqs of length 300 bp)
    num_seqs_in_set = 100
    len_seq = 300

    geom_rvs = geom.rvs(0.75, size=num_seqs_in_set, loc=-1) #sym2=0.75, sym3=.7. Originally 0.5
    
    set1 = MGlib.random_DNA(len_seq,{'A': 0.3,'C': 0.2,'G': 0.2,'T': 0.3},num_seqs_in_set)
    
    # sample large number of sites from motif
    pmot1 = MGlib.sample_motif(mot, num_seqs_in_set) 

    if mode=="positive":
        #insert sites in sequences
        e=0
        while (e<len(set1)):
            # edit sequence to include random site(s)
            # determine number of sites per geometric distribution
            num_sites = geom_rvs[e]
            new_sites = ""
            for j in range(0, num_sites): 
                new_sites += random.choice(pmot1)
            if len(new_sites) > len_seq:
              new_sites = new_sites[:len_seq]
            set1[e] = new_sites + set1[e][len(new_sites):]
            e = e+1

    set2=set1
    return set2


def get_cog_type(which_one, ncutoff):
    ctype = "negative"
    if which_one >= ncutoff:
      ctype = "positive"
    return ctype
    
def compute_log_l(score_set, the_n_g, the_n_m, alpha):
    sumll = 0.0
    for score_list in score_set: # for each of the 100 sequences
        #compute the sum of log likelihood of each score array
        lpd_fs = np.log(alpha*the_n_m.pdf(score_list) + (1-alpha)*the_n_g.pdf(score_list))
        sumll += sum(lpd_fs)
    return sumll

def sample(n, xs, replace=True):
    """Samples n objects from the list xs."""
    if replace:
        return [random.choice(xs) for i in range(n)]
    else:
        ys = list(xs[:])
        samp = []
        for i in range(n):
            y = random.choice(ys)
            samp.append(y)
            ys.remove(y)
        return samp


def permute_pssm(the_pssm):
    arr = the_pssm.pssm
    arr2 = the_pssm.pssm
    nums = [i for i in range(len(arr[0]))]
    pnums = sample(len(nums), nums, replace=False)
    for n in nums:
        for r in range(0,4):
            arr2[r][n] = arr[r][pnums[n]]
    return arr2

def permute_sites(sites):
    """Permutes columns of the binding motif."""
    sites = sites[:]
    l = len(sites[0])

    p = sample(l, range(l), replace=False)
    for i, site in enumerate(sites):
        sites[i] = "".join(site[p[j]] for j in p)

    return sites
    
def permute_motif(cur_motif_sites):
    new_sites = permute_sites(cur_motif_sites)
    new_motif_sites = []
    for new_site in new_sites:
        new_motif_sites.append(Seq(new_site,IUPAC.unambiguous_dna))
    new_motif = motifs.create(new_motif_sites) 
    new_motif.pseudocounts=1
    new_motif.background=None
    return new_motif
    
def sym_permute_sites(sites):
    sites = sites[:]
    # We know CsoR is symmetrical around pos 8, so permute 0..7
    msize = len(sites[0])-1 # Want 16 here
    l = 8
    p_half = sample(l, range(l), replace=False)
    p = []
    for num in p_half:
      p.append(num) # save first half
    p.append(l)
    for num in reversed(p_half):
      p.append(msize-num) # save 2nd half

    for i, site in enumerate(sites):
        sites[i] = "".join(site[p[j]] for j in p)
    return sites
    
def sym_permute_motif(cur_motif_sites):
    new_sites = sym_permute_sites(cur_motif_sites)
    new_motif_sites = []
    for new_site in new_sites:
        new_motif_sites.append(Seq(new_site,IUPAC.unambiguous_dna))
    new_motif = motifs.create(new_motif_sites) 
    new_motif.pseudocounts=1
    new_motif.background=None
    return new_motif
  
def compute_p_val(ll_list, true_ll):
  # Sort the ll_list
  # See where the true_ll falls
  ll_list_sorted = sorted(ll_list, reverse=True)
  pos = len(ll_list)+1
  for i in range(0, len(ll_list_sorted)):
    # Go from largest value to smallest
    val = ll_list_sorted[i]
    if val <= true_ll:
      # True ll would go here in the list
      pos = i+1 # Since i starts counting with 0
      break

  pval = float(pos) / float((len(ll_list)+1))
  #print "PVAL", pval, true_ll, pos, ll_list_sorted
  return pval


def main():
    ###############################################################################
    #set default parameters
    motif_filename="CsoR.txt"   #input file
    out_filename="cog_exp_sym2_c"   #prefix for output
    verbose=0                   #verbose mode
    alpha=1.0/300.0             #mixing ratio for regulated model
    rproms=3.0                  #number of regulated promoters [prior]
    tproms=1811.0               #total number of promoters in genome [prior]

    # control number of cogs and number of permutations
    num_cogs =  10000
    neg_cutoff = 9900 # Cog #'s less than this are negative
    num_perms = 100
    cog_sample_size = 1000

    #verbose
    if verbose: print "Using: ", motif_filename, " as input"
    if verbose: print "Writing to (suffix): ", "[void]" if out_filename==""\
    else out_filename
    
    #open file for ouput
    try:
        out_file = open(out_filename + str(num_cogs)+"_s"+str(cog_sample_size)+"_p"+str(num_perms)+".csv","w")
    except (IOError, OSError) as file_open_exception:
        print "*** Something went wrong while opening the output file"
        print "*** Error: ", file_open_exception.errno, " - ",\
                             file_open_exception.strerror        
        sys.exit()      
        
    #compute priors
    PR=rproms/tproms               #prior probability of regulation
    PB=1.0-PR                      #prior probability of non-regulation
    PPR=PB/PR                      #prior probability ratio
        

    
    # read motif and assign 0.25 pseudocounts to PSWM
    # also assign background uniform distribution for the PSSM (default)
    mot = MGlib.read_motif(motif_filename)
    mot.pseudocounts=1
    mot.background=None

    # save the pssm for the motif and the reverse complement
    #(so that they are not recalculated everytime we invoke motif.pssm)
    pssm = mot.pssm
    rpssm = pssm.reverse_complement()

    
    # Save the motif itself as a list of strings for later permuting
    motif_sites = []
    num_motif_sites = len(mot.instances)
    for i in range(num_motif_sites):
        motif_sites.append(str(mot.instances[i]))
    

    
    random.seed(None)
         
 
    # Create the COGS
    all_cogs = []
    the_neg_seqs = []
    neg_cog_nums = [i for i in range(0, neg_cutoff)]
    ran_neg_cog_nums = sample(cog_sample_size, neg_cog_nums, replace=False)
    cog_file = open("seqs_sym2_"+str(num_cogs)+"_s"+str(cog_sample_size)+"_p"+str(num_perms)+".csv","w")
    for i in range(0,num_cogs):
        label = get_cog_type(i, neg_cutoff)
        #print "Create cog #", i, label
        cur_cog = create_COG(label, mot)
        all_cogs.append(cur_cog)
        if i in ran_neg_cog_nums:
            # A negatively regulated cog
            for s in cur_cog:
                the_neg_seqs.append(s)
                cog_file.write("%d,%s\n" % (i,s))
        else:
            for s in cur_cog:
                cog_file.write("%d,%s\n" % (i,s))
    cog_file.close()
  
   
    # compute softmax scores for sampled background sequences
    gscr = MGlib.esfmax_score_seqs(the_neg_seqs,pssm,rpssm) 
    # compute softmax scores for motif sequences
    mscr = MGlib.esfmax_score_seqs(mot.instances,pssm,rpssm)
        
    # get normal distributions for background and motif
    mean_gscr = mean(gscr)
    std_gscr = std(gscr)
    n_g=norm(mean_gscr, std_gscr)
    mean_mscr = mean(mscr)
    std_mscr = std(mscr)
    n_m=norm(mean(mscr), std(mscr))

    smeans_file = open("smeans_stds_sym2_"+str(num_cogs)+"_s"+str(cog_sample_size)+"_p"+str(num_perms)+".csv", "w")
    smeans_file.write("PSSM n_g,%13.10f,%13.10f\n" % (mean_gscr,std_gscr))
    smeans_file.write("PSSM n_m,%13.10f,%13.10f\n" % (mean_mscr,std_mscr))
    
    # Create the permuted pssm and n_m and n_g for the permutation tests
    new_pssm_list = []
    rnew_pssm_list = []
    n_m_perms = []
    n_g_perms = []

    for j in range(0, num_perms):
        #print "\n***************** Create permutation #", j
        # permute the columns of the motif
        new_mot = sym_permute_motif(motif_sites)
        new_pssm = new_mot.pssm #
        rnew_pssm = new_pssm.reverse_complement()
        new_pssm_list.append(new_pssm)
        rnew_pssm_list.append(rnew_pssm)
        # compute score for the negative sequences
        gscr = MGlib.esfmax_score_seqs(the_neg_seqs,new_pssm,rnew_pssm)
        mean_gscr = mean(gscr)
        std_gscr = std(gscr)
        # compute softmax scores for new motif sequences
        mscr = MGlib.esfmax_score_seqs(new_mot.instances,new_pssm,rnew_pssm)
        mean_mscr = mean(mscr)
        std_mscr = std(mscr)
        smeans_file.write("PermPSSM n_g,%13.10f,%13.10f\n" % (mean_gscr,std_gscr))
        smeans_file.write("PermPSSM n_m,%13.10f,%13.10f\n" % (mean_mscr,std_mscr))

        # get normal distributions for background and motif
        n_g_temp=norm(mean_gscr, std_gscr)
        n_g_perms.append(n_g_temp)
        n_m_temp=norm(mean_mscr, std_mscr)
        n_m_perms.append(n_m_temp)
    smeans_file.close()
    
    # write csv header
    out_file.write('COG Num,Pos/Neg Regulated,Posterior,LogLikelihood,True Model LL,LL Pval\n')
      
    # For each cog, do the posterior calculation and the permutation tests
    for i in range(0,num_cogs):
        label = get_cog_type(i, neg_cutoff)
        #print "Test Cog:", i,label
        # The original posterior computation
        #compute softmax scores for sequences in dataset
        scrs=MGlib.esfmax_score_seqs(all_cogs[i],pssm,rpssm)
        #print np.min(scrs[0]), np.max(scrs[0])
        # Compute posterior 
        # get log-likelihoods for sequences in dataset        
        llrs=MGlib.ll_ratios(scrs,n_g,n_m,alpha)
        # get per-sequence posterior for the sequences in dataset
        fpost=MGlib.PostP(llrs,PPR,0)
    
        true_model_ll = compute_log_l(scrs, n_g, n_m, alpha) 
        
        
        #####################################
        # Permutation test
        log_ls = []
        for j in range(0, num_perms):
            #print " ... perm test", j
            # Compute score and log likelihood for each permutation.
            scrs=MGlib.esfmax_score_seqs(all_cogs[i],new_pssm_list[j],rnew_pssm_list[j])
            log_l = compute_log_l(scrs, n_g_perms[j], n_m_perms[j], alpha)
            log_ls.append(log_l)
    
        rev_pval = compute_p_val(log_ls, true_model_ll)
        pval = 1.0 - rev_pval
        out_file.write("%d,%s,%10.7f,%10.7f,%10.7f,%10.7f\n" % (i, label, fpost, rev_pval, true_model_ll, pval))

     
    out_file.close()

#if __name__ == "__main__": 
main()
    
