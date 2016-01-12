# Get the information for a particular TF and COG from the "*_genes.csv" files.
# Theta of interest is hardcoded below.
# NOTE: will reverse complement the - strand promoters

import sys
import os
import glob

def get_sample_from_fname(fname):
  basename = os.path.basename(fname)
  items = basename.split("_")
  spos = 0 
  return items[spos]
  
the_dir = sys.argv[1]
TF = sys.argv[2]
cog = sys.argv[3]
theta = "-6.000000" #"-2.000000"

file_to_glob = os.path.join(the_dir, "*_genes.csv")
#print file_to_glob
file_array = glob.glob(file_to_glob)

cpos = 0 # cog or nog
thpos = 1 # theta value
gpos = 3  # gene name
ppos = 4 # promoter sequence (not reverse complemented for - strands yet)
stpos = 5 # strand (+ or -)
npos = 6   # count

# http://stackoverflow.com/questions/19570800/reverse-complement-dna
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a'}[B] for B in x][::-1])
for f in file_array:
  sname = get_sample_from_fname(f)
  #print sname
  with open(f, "r") as myfile:
    inlines = myfile.readlines()
  #print sname, len(inlines)
  for l in inlines:
    # No header line in these
    line = l.strip()
    if len(line) < 1:
      continue
    items = line.split(',')
    if items[cpos] == cog and items[thpos] == theta:
      if items[stpos] == "+":
        prom_seq = items[ppos]
      else:
        prom_seq = revcompl(items[ppos])
        
      out_str = sname+","+items[gpos]+","+TF+","+prom_seq+","+items[stpos]+","+items[npos]
      print out_str