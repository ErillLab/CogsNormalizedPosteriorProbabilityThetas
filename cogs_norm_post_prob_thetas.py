# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:32:51 2015

NOTE: this file filters by taxonomy -- Firmicutes or Gammaproteobacteria
-- and by head completeness (Complete or Lack 3'-end)
@author: Based on Talmo's original analysis code
"""

import sys
import os
#sys.path.append("../..")
from igc_pipeline import *
import pandas as pd
import numpy as np
import scipy
import scipy.stats

# Configuration
TF = sys.argv[1]

firmi = True # change to gamma if want LexA gamma
skip=True # to skip short promoters
out_extra = ""

# test_cogs identifies cogs for debug print
# out_cogs identifies cogs to be written out to a file
num_reg_cogs = 3.0 # Use 3 for CsoR and YrkD
total_num_cogs = 1811.0 # Use 1811 for CsoR and YrkD
if TF == "YrkD":
  bsite_fname = "Bacillales_YrkD.txt"
  bg_mu_val=-17.186085
  bg_sigma_val=7.213475
  test_cogs = ["COG1937", "COG2217", "NOG216437", "COG3550", "NOG36059", "NOG72983", "NOG63000", "NOG02494"] #["COG1937"] # YrkD
  out_cogs = ["COG1937", "COG2217", "NOG216437", "COG3550", "NOG36059", "NOG72983", "NOG63000", "NOG02494"] 
elif TF=="LexA":
  num_reg_cogs = 30.0 # Use 3 for CsoR and YrkD
  total_num_cogs = 2000.0 # Use 1811 for CsoR and YrkD
  if firmi:
    bsite_fname = "Firmicutes_LexA.txt" #-17.681878, bg_sigma=8.2267625
    bg_mu_val= -17.681878
    bg_sigma_val=8.2267625
    print "Handling Firmicutes"
    test_cogs = ["COG4277","COG5108", "COG1974", "NOG38629", "NOG79879", "NOG46871", "NOG101678"] #["COG0799","NOG46871","COG5108", "COG1974"]
    out_cogs = ["COG4277","COG5108", "COG1974", "NOG38629", "NOG79879", "NOG46871", "NOG101678"] #["COG0799","NOG46871","COG5108", "COG1974"]
  else:
    bsite_fname = "GammaProteobacteria_LexA.txt" #-21.660313, bg_sigma=8.476820
    bg_mu_val=-21.660313
    bg_sigma_val=8.476820
    print "Handling Gammaproteobacteria"
    out_cogs = ["COG1974", "COG3070","NOG82272"] #["COG1546", "COG1974", "COG0799","NOG46871","COG5108"]
    test_cogs = ["COG1974", "COG3070","NOG82272"] #["COG1546", "COG1974", "COG0799","NOG46871","COG5108"]
  #test_cogs = ["COG1974"] # LexA
  #out_cogs = ["COG1974"]
elif TF=="CsoR":
  bsite_fname = "Bacillales_CsoR.txt"
  bg_mu_val=-16.454468
  bg_sigma_val=6.663047
  test_cogs = ["COG1937", "COG2217", "NOG218972", "NOG72602", "NOG109008", "COG2836", "COG4633", "NOG81268", "NOG216622", "NOG10047", "NOG204883", "NOG89021", "NOG216557"] #["COG1937"] # CsoR 
  out_cogs = ["COG1937", "COG2217", "NOG218972", "NOG72602", "NOG109008", "COG2836", "COG4633", "NOG81268", "NOG216622", "NOG10047", "NOG204883", "NOG89021", "NOG216557"] 
elif TF=="CopY":
  bsite_fname = "CopY.txt"
  bg_mu_val=-27.392480
  bg_sigma_val=10.159172
  test_cogs = ["COG0745", "COG3682", "COG3142", "COG2132", "COG2372", "COG2608", "COG2847", "COG3667", "COG1276", "COG1277", "COG4454"] # 
  out_cogs = ["COG0745", "COG3682", "COG3142", "COG2132", "COG2372", "COG2608", "COG2847", "COG3667", "COG1276", "COG1277", "COG4454"] 
else:
  print "Need to set up this TF first"
  sys.exit(0)

if not firmi:
  out_extra = "_gamma"


change_score_path("4TB")
    
pssm = PSSMScorer(binding_sites_path + bsite_fname, TF)
pssm.initialize_estimator(bg_mu=bg_mu_val, bg_sigma=bg_sigma_val)

pname = pssm.name

self_scores = pssm.score_self(True)
mean_self = np.mean(self_scores)
mean_std = np.std(self_scores)

prob_f = num_reg_cogs/total_num_cogs
prob_b = 1.0 - prob_f

samples = get_all_with_scores(pssm)

# samples = ["MH0012"] # for testing

min_plen = 60 - pssm.length
if skip:
  skipstr="_skipshort"
  print "Skip short promoters", min_plen
else:
  skipstr = ""
  print "Not skipping short promoters"

# Specify the output directory for the run
out_dir = "/home/cuda/2TB/Output/"+TF+"_cog_norm"+out_extra+skipstr

if not os.path.exists(out_dir):
  os.mkdir(out_dir)

print "Mean self scores:", mean_self, "stddev:", mean_std
cdf_background_dist = scipy.stats.distributions.norm(bg_mu_val, bg_sigma_val).cdf
alpha = 1.0/300.0 # Note same as in PSSMScorer
cdf_foreground_dist = scipy.stats.distributions.norm(mean_self, mean_std).cdf
sdev_mults = [-0.5, -1.0, -1.5, -2.0, -3.0, -4.0, -6.0, -8.0]
#sdev_mults = [-2.0]
thetas = {}
use_prom_len = 300.0
for mult_val in sdev_mults:
    theta = mean_self + (mean_std * mult_val) # use neg mult val to go below mean
    print "Multiplier:", mult_val, "Theta:", theta
    cdf_background_val = cdf_background_dist(theta)
    print "CDF background value for theta:", cdf_background_val
    cdf_foreground_val = alpha * cdf_foreground_dist(theta) + (1 - alpha) * cdf_background_dist(theta)
    print "CDF foreground value for theta:", cdf_foreground_val
    U = np.power(cdf_background_val, use_prom_len)
    OneMinusU = 1.0 - U
    Up = np.power(cdf_foreground_val, use_prom_len)
    OneMinusUp = 1.0 - Up
    norm_p_f = (num_reg_cogs * OneMinusUp)/((num_reg_cogs*OneMinusUp) + (total_num_cogs*OneMinusU))
    norm_p_b = 1.0 - norm_p_f
    print "U:", U, "1-U", OneMinusU, "U'", Up, "1-U'", OneMinusUp
    print "NormP(f)",norm_p_f, "normP(B)", norm_p_b
    prior = norm_p_b/norm_p_f
    thetas[mult_val] = {"theta":theta, "cdf_bg":cdf_background_val, "cdf_fg":cdf_foreground_val, "prior":prior}
    print thetas[mult_val]

scnt = 0
for sample in samples:
  scnt += 1
  if not sample[0:2]=="MH": #    CHANGE WHICH MH or non-MH Samples here  !!!!!!!!!!!!!!!
      continue
  print "sample:", sample
  try:
    # Load data
    genes = load_sample_genes(sample)
    operons = get_operons(sample)

    # genes2operon is a series that maps gene names to operons
    # genes2operon has genename and a number like 0, 1, 2...
    genes2operon = get_genes2operon(genes, operons)

    scores = get_sample_scores(sample, pssm)
    LL_ratios = pd.Series()

    # Filter by head completeness. Need to remove both Lack 5'-end and Lack both ends
    # Having trouble with compound if test for Series (which is what operons.genes is)
    ol5 = operons[operons.head_completeness!="Lack 5'-end"]
    ol5e = ol5[ol5.head_completeness!="Lack both ends"]
    #genes = genes.loc[np.unique(np.hstack(operons.genes[operons.head_completeness != "Lack 5'-end"]))].dropna(how="all")
    genes2 = genes.loc[np.unique(np.hstack(ol5e.genes))].dropna(how="all")
    genes = genes2.drop_duplicates()
    genes.index.name = "gene_name"
    # Filter by taxonomy
    if firmi:
      genes = genes[(genes.phylum == "Bacteroidetes") ]#| (genes.phylum == "Actinobacteria")] [(genes.phylum == "Firmicutes") ]
    else:
      genes = genes[genes["class"] == "Gammaproteobacteria"]
    #tax_tf_genes = sum(genes.eggNOG == the_cog)
  
    # Group genes by COG
    grouped = genes.reset_index().groupby("eggNOG")

    # Talmo's original code for computing the posterior probability is nicely vectorized.
    # However, there was a problem with not finding known regulators properly.
    # In the process of debugging this problem, the vectorized code was unwound into
    # the loops below. Also, individual values were output into CSV files to eaable performing
    # separate analyses in Excel or other for these values.
    # At the same time, it was decided to use cutoff values for excluding some hits.
    # Ultimately, the bug was not in this analysis code but in how operons were scored.
    # However, by then everything was working so the code was left as loops instead of re-vectorizing it.

    # Compute log-likelihood ratios and posteriors - Talmo's original code, for reference.
    # This is how the code originally looked.
    #get_LL_ratio = lambda gene_names: pssm.LL_ratio(np.hstack(scores.loc[genes2operon[gene_names]]))
    #sample_ratios = grouped["gene_name"].agg(get_LL_ratio)
    #idx = LL_ratios.index.union(sample_ratios.index)
    #LL_ratios = LL_ratios.get(idx).fillna(1) * sample_ratios.get(idx).fillna(1)
    #cog_probs = LL_ratios.map(lambda LLR: 1 / (1 + LLR * pssm.pb / pssm.pf))

    # Files to contain the normalized cog probabilities per sample and
    # the genes information. This data will be used in the post-processing to determine posteriors for
    # individual cogs and to look at specific genes..
    out_fname = os.path.join(out_dir, sample+"_"+TF+"_norm_cog_probs.csv")
    ofp = open(out_fname, "w")
    out_gname = os.path.join(out_dir, sample+"_"+TF+"_genes.csv")
    ogp = open(out_gname, "w")
    # Process each cog
    for the_cog, the_genes in grouped:
      if the_cog in test_cogs:
        print "COG:", the_cog
      #the_genes = grouped.get_group(the_cog)
      the_gene_names = the_genes["gene_name"] # a Series

      for mult_val in sdev_mults:
          # Since we are decrementing total_genes if skipping short
          # need to reset its value for each multiplier
          total_genes = len(the_gene_names)
          # Set the normalization tally variables to 0
          num_mapping_promoters = 0 # |PC|      
          total_log_OneMinusU = 0.0
          total_log_OneMinusUp = 0.0
          total_log_PB_PF = 0.0
          normalized_ratio = 1.0 # In case no good promoters
          cpprob = 1/(1+normalized_ratio*thetas[mult_val]["prior"]) # default
          # Output information for each gene for this cog individually
          for a_gene in the_gene_names:
            # gene2operon is a series returned from get_genes2operon
            # Can reference items in a series by [pos] and ['name']. Also [pos,pos,pos].
            g_op = genes2operon[a_gene] # returns the operon for the a_gene
            gene_scores = scores.loc[g_op] # returning array of scores from the sliding window over promoter
            len_Pc = len(gene_scores) # Number of scores in the promoter is basically its length
            if skip and (len_Pc <= min_plen):
              total_genes = total_genes - 1 # since ignoring short, don't count toward total
              if the_cog in test_cogs:
                print "Skip gene", a_gene, len_Pc
              continue
            promseq = operons[g_op:g_op+1]["promoter_seq"]
            the_prom = promseq[g_op]
            strandseq = operons[g_op:g_op+1]["strand"]
            the_strand = strandseq[g_op]
            good_scores = gene_scores[gene_scores >= thetas[mult_val]["theta"]]
            len_Pc_good = len(good_scores)
            max_score = max(gene_scores)
            mspos_str = str(np.argmax(gene_scores))+",0"
            ms_str = str(max_score)
            glr = 1.0
            gcpO = 1.0/(1.0 + glr*pssm.pb/pssm.pf)
            gcpP = 1.0/(1.0 + glr*thetas[mult_val]["prior"])
            if the_cog in test_cogs:
              print "  Gene:", a_gene
              print "    Num scores (approx promoter len):", len_Pc
              print "    Num good scores:", len_Pc_good
              #print "    Promoter:", the_prom, the_strand
            
            if len_Pc_good > 0:
              if len_Pc_good > 1:
                ms_str = ""
                for ascore in good_scores:
                  ms_str += str(ascore)+";"
                if len(ms_str) > 1:
                  ms_str = ms_str[:-1]
                winfo = np.where(gene_scores >= thetas[mult_val]["theta"])
                the_pos = winfo[0]
                if len(the_pos) > 1:
                  mspos_str = ""
                  for posval in the_pos:
                    mspos_str += str(posval)+","
                  if len(mspos_str) > 1:
                    mspos_str = mspos_str[:-1]
              # There are scores >= theta, so work with this promoter
              num_mapping_promoters += 1
    
              # The following code lines follow the normalized probability discussion.
              # Values are computed following the original equations for comparison to the normalized.
              U = np.power(thetas[mult_val]["cdf_bg"], len_Pc)
              OneMinusU = 1.0 - U
    
              Up = np.power(thetas[mult_val]["cdf_fg"], len_Pc)
              OneMinusUp = 1.0 - Up
              total_log_OneMinusU += np.log(OneMinusU) # This for the (1-U)^^|PC|
              total_log_OneMinusUp += np.log(OneMinusUp) # This for the (1-U')^^|PC|
              if the_cog in test_cogs:
                print "    U :", U, "1-U :", OneMinusU
                print "    U':", Up, "1-U':", OneMinusUp
                print "    Total log(1-U):", total_log_OneMinusU, "Total log(1-U'):", total_log_OneMinusUp
              # LL_ratio equation: np.exp(np.sum(np.log(self.L_b(scores)) - np.log(self.L_f(scores))))        
              # Reconstruct 
              bscores = pssm.L_b(gene_scores) # P(Pc|B) -- use all scores in promoter
              fscores = pssm.L_f(gene_scores) # P(Pc|F)
              logbscores = np.log(bscores)
              logfscores = np.log(fscores)
              ldiff = logbscores - logfscores # This is the division P(Pc|B)/P(Pc|F), may be an array
              ldiff_llratio = np.exp(np.sum(ldiff)) # Original equation for LLR
              glr = ldiff_llratio
              ldiff_cogprobO = 1/(1+ldiff_llratio*pssm.pb/pssm.pf)
              ldiff_cogprobP = 1/(1+ldiff_llratio*thetas[mult_val]["prior"])
              gcpO = ldiff_cogprobO
              gcpP = ldiff_cogprobP
              log_ratio = np.sum(ldiff)
              total_log_PB_PF += log_ratio # For the multiplication of the P(Pc|B)/P(Pc|F) ratio
              if the_cog in test_cogs:
                # prior used to be pssm.pb/pssm.pf
                # This is on a per-gene basis, not per-cog, so not using
                norm_p_f = (num_reg_cogs * OneMinusUp)/((num_reg_cogs*OneMinusUp) + (total_num_cogs*OneMinusU))
                norm_p_b = 1.0 - norm_p_f
                ldiff_cogprob = 1/(1+ldiff_llratio*norm_p_b/norm_p_f)

                #print "    Log backgr scores:", logbscores
                #print "    Log foregr scores:", logfscores
                #print "    Log difference b-f:", ldiff
                print "    Multiplier", mult_val, "Theta:", thetas[mult_val]["theta"], "Prior Rat", thetas[mult_val]["prior"]
                print "    Norm p_f:", norm_p_f, "norm_p_b:", norm_p_b
                print "    Sum of log diff (log ratio):", log_ratio
                print "    Total log ratio thus far:", total_log_PB_PF
                print "    Original LLR for gene:", ldiff_llratio
                print "    Original cog post prob:", ldiff_cogprobO
                print "    Cog post prob (gene prior ratio):", ldiff_cogprob
                print "    Cog post prob (theta prior ratio):", ldiff_cogprobP
                val = ldiff_llratio * (OneMinusUp/OneMinusU)
                print "    Alt val:", val
                #print "    Good scores:", good_scores
            if the_cog in out_cogs:
              ogp.write("%s,%f,%12.10f,%s,%s,%s,%d,%s,%12.10f,%12.10f,%12.10f,%s\n" % (the_cog, mult_val, thetas[mult_val]["theta"], a_gene, the_prom, the_strand, len_Pc_good, ms_str, glr, gcpO, gcpP, mspos_str))
          # end loop genes for the cog
          # Handle this cog
          if num_mapping_promoters > 0: # we have something that is above the threshold.
            cog_OneMinusU = np.exp(total_log_OneMinusU)
            cog_OneMinusUp = np.exp(total_log_OneMinusUp)
            cog_PB_PF = np.exp(total_log_PB_PF)
            normalized_ratio = cog_PB_PF * (cog_OneMinusUp / cog_OneMinusU)
            cpprob = 1/(1+normalized_ratio*thetas[mult_val]["prior"])
            if the_cog in test_cogs:
              print the_cog, "# mapping promoters:", num_mapping_promoters
              print "Cog's OneMinusU':", cog_OneMinusUp
              print "Cog's OneMinusU :", cog_OneMinusU
              print "Cog LB/LF:", cog_PB_PF
              print "Normalized ratio:", normalized_ratio
              print "Prior ratio based on theta:", thetas[mult_val]["prior"]
              print "Cog's normalized post prob:", cpprob
              print " " 
          ofp.write("%s,%f,%12.10f,%12.10f,%12.10f,%12.10f,%d,%d\n" % (the_cog, mult_val, thetas[mult_val]["theta"],thetas[mult_val]["prior"],normalized_ratio, cpprob,total_genes,num_mapping_promoters))
      # end theta multipliers loop
    # end cog loop
    ofp.close()
    ogp.close()
  except:
    print "Error:", sys.exc_info()[1]
