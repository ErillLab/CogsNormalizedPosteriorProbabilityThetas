import os
import sys
import time

gfiles_dir =  "../csor/input/CsoR_cog_norm_skipshort" 
out_dir =  "../csor/output/csor_prom_seqsrun1" 
TF =  "CsoR" 
val_type = "-6.0" # "-2.0" # match to gather_gene_info.py.
type = "firmi" 
cog = sys.argv[1] 

init_fname = TF+"_"+type+"_genes_"+cog+"_theta"+val_type+".csv"
out_fname1 = os.path.join(out_dir, init_fname)
sys_cmd = "python gather_gene_info.py "+gfiles_dir+" "+TF+" "+cog+" > "+out_fname1
print sys_cmd
os.system(sys_cmd)

time.sleep(2)

second_fname = "fastas_"+init_fname
out_fname2 = os.path.join(out_dir, second_fname)
sys_cmd = "python get_gene_fastas_alt.py " + out_fname1
print sys_cmd
os.system(sys_cmd)

time.sleep(2)
# Need to rename the output file
print "rename", second_fname, "to", out_fname2
os.rename(second_fname, out_fname2)

time.sleep(2)
out_fname3 = os.path.join(out_dir, TF+"_"+type+"_"+cog+"_gene_proms_info_theta"+val_type+".csv")
sys_cmd = "python obtain_prom_seqs.py " + out_fname2+" "+out_fname1+" > "+out_fname3
print sys_cmd
os.system(sys_cmd)
