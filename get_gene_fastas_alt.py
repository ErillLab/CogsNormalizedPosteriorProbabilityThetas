import os
import sys
import glob
# Go through the CSV file created with ___
# It has samples, genes (one per line), the op #, op position, cog prob

#base_path = "/UMBC/PhD/ErillRotation/metagenomics/"
base_path = "E:/ErillRotation/metagenomics/"
IGC_path = base_path + "IGC/"
orfs_fasta_path = IGC_path + "5.IndividualORFs/"
o_suffix = ".fna" # new style

fname = sys.argv[1]
base_name = os.path.basename(fname)
lcnt = 0
lstart = 1 # NO header
with open(fname, "r") as myfile:
  inlines = myfile.readlines() # file's not that big
  
spos = 0
gpos = 1
samples = []
sample_genes = {}

for l in inlines:
  lcnt += 1
  if lcnt < lstart:
    continue
  line = l.strip()
  if len(line) < 1:
    continue
  items = line.split(",")
  sample_name = items[spos]
  gene = items[gpos]
  #print sample_name, gene
  if sample_name not in samples:
    samples.append(sample_name)
  if sample_name not in sample_genes:
    #print "Add new sg"
    sample_genes[sample_name] = [gene]
  else:
    curlist = sample_genes[sample_name]
    #print "cur list:", curlist
    curlist.append(gene)
    sample_genes[sample_name] = curlist
  #print "sample_genes is", sample_genes[sample_name]

ofp = open("fastas_"+base_name, "w")
for sample_name in samples:
  fasta_file = os.path.join(orfs_fasta_path, sample_name+o_suffix)
  print "Open", fasta_file
  if not os.path.exists(fasta_file):
    print " Fasta not found"
    continue
  genes = sample_genes[sample_name] # Which genes we care about
  genes_printed = []
  out_str = ""
  with open(fasta_file, "r") as myfile:
    do_print = False
    for il in myfile:
      inline = il.strip() # term or phrase
      if len(inline) < 1:
        continue

      if inline[0] == ">":
        if len(out_str) > 0:
          ofp.write("%s\n" % out_str)
        out_str = ""
        # Some fasta files have a tab first
        spos = inline.find("\t")
        if spos < 0:
          spos = inline.find(" ")
        if spos > -1:
          gname = inline[1:spos]
          if gname in genes:
            out_str = inline+"\n"
            do_print = True
            genes_printed.append(gname)
          else:
            do_print = False

      elif do_print:
        out_str += inline
        #do_print = False # The Stool fastas are multi line
  if len(out_str) > 0:
    ofp.write("%s\n" % out_str) # last one  
  # Check all genes picked up
  for g in genes:
    if g not in genes_printed:
      print sample_name, g, "not found"
