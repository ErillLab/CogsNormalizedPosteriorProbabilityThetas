import os
import sys

fastas = {}

def read_in_fasta_file(fasta_file):
  with open(fasta_file, "r") as myfile:
    do_print = False
    out_str = ""
    for il in myfile:
      inline = il.strip() # term or phrase
      if len(inline) < 1:
        continue
      #print "Line:", inline
      if inline[0] == ">":
        #print "> out_str?", out_str, len(out_str)
        if len(out_str) > 0:
          #print gname, "SAVE", out_str
          fastas[gname]= out_str

        out_str = ""
        # Some fasta files have a tab first
        spos = inline.find("\t")
        if spos < 0:
          spos = inline.find(" ")
        if spos > -1:
          gname = inline[1:spos]
          #print "GN", gname
          do_print = True

      elif do_print:
        out_str += inline
        #print "Out:", out_str
        #do_print = False # The Stool fastas are multi line
  if len(out_str) > 0:
    #print "last",gname,"is",out_str
    fastas[gname] = out_str # last one  

fname = sys.argv[1] # Name of file with the fastas
gname1 = sys.argv[2] # Name of file with gene info for YrkD
gname2 = ""
if len(sys.argv) > 3:
  gname2 = sys.argv[3] 

read_in_fasta_file(fname)
#print fastas

# Since gname1 and gname 2 have the same genes, we only need the
# 2nd one to get whether or not it matched above theta or not.

spos = 0
gpos = 1
tpos = 2
ppos = 3
stpos = 4
npos = 5

genes = {}

with open(gname1, "r") as myfile:
  in_lines = myfile.readlines()
  
# Read in the first one and set up the structure

for l in in_lines:
  line = l.strip()
  if len(line) < 1:
    continue
  items = line.split(",")
  gname = items[gpos]
  if gname not in fastas:
    print "ERROR: no fasta for", gname
    sys.exit(0)
  else:
    num = int(items[npos])
    if num > 0:
      TF = items[tpos]
    else:
      TF = "none" # may be overwritten by 2nd file
    genes[gname] = [TF,items[stpos],items[ppos]]
	
if len(gname2) > 0:
  # 2nd file
  with open(gname2, "r") as myfile:
    in_lines = myfile.readlines()

  for l in in_lines:
    line = l.strip()
    if len(line) < 1:
      continue
    items = line.split(",")
    gname = items[gpos]
    if gname not in genes:
      print "ERROR: new gene", gname
      sys.exit(0)
    else:
      num = int(items[npos])
      if num > 0:
        TF = items[tpos]
        # Overwrite None with new TF
        genes[gname] = [TF,items[stpos],items[ppos]]
	  
gkeys = genes.keys()
gkeys.sort()
out_str = "Gene,TF,Strand,Promoter+Gene beginning,Promoter,Gene beginning,Original Promoter,Original Gene"
print out_str
for gkey in gkeys:
  ginfo = genes[gkey]
  out_str = gkey+","+ ginfo[0]+","+ginfo[1]+"," + ginfo[2].upper()+fastas[gkey][0:50].upper()+ ","+ginfo[2].upper()+","+fastas[gkey][0:50].upper()+","+ginfo[2]+","+fastas[gkey]
  print out_str
