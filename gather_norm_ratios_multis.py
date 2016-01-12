# This script processes the files matching *norm_cog_probs.csv in an input directory
# and extracts each COG's information separately for each theta.
# The normalized cog probabilities are logged and summed. The normalized prior is obtained from the CSV.
# For the COG for a given theta across all samples, the sum of the log ratio, the normalized ratio,
# the posterior probability, the total genes, and the total genes passing the theta.

import os
import sys
import glob
import numpy as np
#import math as np
from operator import itemgetter

def get_sample_from_fname(fname):
  basename = os.path.basename(fname)
  items = basename.split("_")
  spos = 0 #1
  return items[spos]

def compute_cog_prob(sum_ln, theta_prior_rat):
  llr = np.exp(sum_ln)
  cp = 1/(1+llr * theta_prior_rat)
  return cp

  
###############################################
in_dir = sys.argv[1] # Where the files are
nroot = sys.argv[2] # "root" for output file name
TF = sys.argv[3]
the_cog = sys.argv[4]
mh_only = True
theta_prior_rats = {}

ftofind = "*norm_cog_probs.csv"
if mh_only:
  ftofind ="MH"+ftofind
file_to_glob = os.path.join(in_dir, ftofind)
file_array = glob.glob(file_to_glob)

theta_llrs = {} # by theta (well, by the multiplier), then cog, the log of the ratio, summed
theta_sample_counts = {} # by cog
theta_tgcnts = {} 
theta_fgcnts = {}
mult_vals = []
c5 = 1.0 # This was just to test the summing and to compare with multiplying
c_interest = 1.0 # Also just to test the summing and to compare with multipyling

cog_pos = 0
mult_pos = 1  # what multiplier * stdev forcomputing theta
prior_pos = 3 # normalized prior
llr_pos = 4 # Normalized ratio
# pos 5 is the normalized cog posterior probability
tg_pos = 6 # total # genes for cog
fg_pos = 7 # Number genes above theta

for f in file_array:
  sname = get_sample_from_fname(f)
  #print sname

  with open(f, "r") as myfile:
    # Read in the llr data for the sample
    inlines = myfile.readlines()
  print sname, len(inlines)
  for l in inlines:
    # No header line in these
    line = l.strip()
    if len(line) < 1:
      continue
    #print line
    items = line.split(',')
    cog = items[cog_pos]
    mult_val = items[mult_pos] # key into 
    if mult_val not in mult_vals:
      mult_vals.append(mult_val)
    llr = float(items[llr_pos])
    tgcnt = int(items[tg_pos])
    fgcnt = int(items[fg_pos])
    if mult_val not in theta_prior_rats:
      theta_prior_rats[mult_val] = float(items[prior_pos])
      print "Adding", mult_val, items[prior_pos]
    # We have a few llr that are 0
    if llr <= 0.0000000000000001:
      llr_log = -36.0
      print "RESET llr_log:", sname, cog, llr
    else:
      #print "Do np.log next of", llr
      llr_log = np.log(llr)
    #print sname, "cog:", cog, "mult_val:", mult_val, "llr:", llr
    if mult_val not in theta_llrs:
      # Create a new llrs
      llrs = {}
      #print "  new llrs"
    else:
      llrs = theta_llrs[mult_val]
      #print "  cur llrs:", llrs
    # For this mult_val, update the cog
    if cog not in llrs:
      llrs[cog] = llr_log
    else:
      llrs[cog] = llrs[cog] + llr_log
    theta_llrs[mult_val] = llrs
    #print "  updated llrs", llrs
    #print "THETA llrs", theta_llrs
 
    if mult_val not in theta_sample_counts:
      sample_counts = {}
    else:
      sample_counts = theta_sample_counts[mult_val]    
    if cog not in sample_counts:
      sample_counts[cog] = 1
    else:
      sample_counts[cog] = sample_counts[cog] + 1
    theta_sample_counts[mult_val] = sample_counts
    #print "TH scnts", theta_sample_counts
    if mult_val not in theta_tgcnts:
      tgcnts = {}
    else:
      tgcnts = theta_tgcnts[mult_val]
    if cog not in tgcnts:
      tgcnts[cog] = tgcnt
    else:
      tgcnts[cog] = tgcnts[cog] + tgcnt
    theta_tgcnts[mult_val] = tgcnts
    #print "TH gcnts:", theta_tgcnts
    if mult_val not in theta_fgcnts:
      fgcnts = {}
    else:
      fgcnts = theta_fgcnts[mult_val]
    if cog not in fgcnts:
      fgcnts[cog] = fgcnt
    else:
      fgcnts[cog] = fgcnts[cog] + fgcnt
    theta_fgcnts[mult_val] = fgcnts
    #print "TH fgcnts", theta_fgcnts
    
    #if cog=="COG0005":
    #  print "--- COG0005: ", cog, llr, llr_log, gcnts[cog]
    #  c5 = c5 * llr
    #  print "---  totals:", llrs[cog], sample_counts[cog], gene_counts[cog], "c5:", c5
    if cog==the_cog:
      print "***",the_cog+": ", cog, llr, llr_log, tgcnt, fgcnt
      c_interest = c_interest * llr
      #print "***  totals:", llrs[cog], sample_counts[cog], gene_counts[cog], the_cog, "val:", c_interest
      print "***  totals:", llrs[cog], sample_counts[cog], the_cog, c_interest, tgcnts[cog], fgcnts[cog]

print "Ready to output now"
for mult_val in mult_vals:
  print "Do multiplier", mult_val
  llrs = theta_llrs[mult_val]
  sample_counts = theta_sample_counts[mult_val]
  tgcnts = theta_tgcnts[mult_val]
  fgcnts = theta_fgcnts[mult_val]
  # Output the cogs. Sort them alphabetically first.
  cog_list = sorted(llrs.keys())
  ofroot = nroot+"_"+TF+"_all_cogs_multis"+mult_val[:4]
  ofname = ofroot+"_summary_norms.csv"
  ofp = open(ofname, "w")
  ofp.write("COG,sum(ln(LLR)),Total Normalized Ratio,Total cog prob new,Total cog prob orig,Total #samples,Total #Genes,Total #Filtered Genes\n")
  for cog in cog_list:
    sum_ln_llrs = llrs[cog]
    cp = 0.0
    cpold = 0.0
    if sum_ln_llrs < 700:
      cp = compute_cog_prob(sum_ln_llrs, theta_prior_rats[mult_val])
      cpold = compute_cog_prob(sum_ln_llrs, 0.99/0.01)
    norm_ratio = 1.0
    try:
      if llrs[cog] > 100:
        # Overflow
        norm_ratio = 2.6e48
      else:
        norm_ratio = np.exp(llrs[cog])
    except:
      print "Exception", llrs[cog]

    outstr = cog+","+str(llrs[cog])+","+str(norm_ratio)+","+str(cp)+","+str(cpold)+","+str(sample_counts[cog])+","+str(tgcnts[cog])+","+str(fgcnts[cog])
    ofp.write("%s\n" % outstr)
  ofp.close()

