# -*- coding: utf-8 -*-
"""
This script handles integrating the results from BLASTing IGC ORFs and/or genes
against the eggNOG database to get COG and taxonomy assignments.

These values are stored in the genes tables for each sample. See 
save_sample_genes() and load_sample_genes() in igc_pipeline for gene table
creation and usage.

Created on Thu Apr 23 18:26:07 2015

@author: Talmo
"""

from igc_pipeline import *
from Bio import Entrez
from tqdm import tqdm
from glob import glob
import os

#%% Configuration
# Paths
eggnog_path = "/home/cuda/2TB/eggNOG/"
members_path = eggnog_path + "NOG.members.txt"
funccats_path = eggnog_path + "NOG.funccat.txt"
blast_out_path = IGC_path + "BLAST/has_cog-eggNOG/"

# NCBI
Entrez.email = "erill@umbc.edu"

#%% Cache data files
t = time()
# Load eggNOG files
eggnog_members = pd.read_table(members_path, index_col="protein name")
eggnog_funccat = pd.read_table(funccats_path, names=["cog", "funccat"], index_col="cog")

# Load taxonomy lookup table
eggnog_tax = pd.read_csv(eggnog_tax_path, index_col="tax_id")

log("Loaded eggNOG data files.", t)

#%% Find not BLASTed
samples = get_all_samples()
blasted = [p.replace(blast_out_path, "").replace(".tbl", "") for p in glob(blast_out_path + "*.tbl") if os.stat(p).st_size > 0]
not_blasted = [sample for sample in samples if sample not in blasted]

#%% Processing functions
def extract_ORFs_with_COGs(sample):
    """ Extracts ORFs for genes that have COG assignment and saves their
    FASTA sequences. """
    sample = get_unique_sample(sample)
    sample_paths = get_sample_paths(sample)
    
    # Load sample data
    genes = load_sample_genes(sample)
    
    # Get ORFs with COGs
    t = time()
    has_cog = genes[genes.eggNOG != "unknown"].index.tolist()
    ORF_seqs = SeqIO.index(sample_paths["ORFs_fasta"], "fasta")
    valid_orfs = [ORF_seqs[orf] for orf in has_cog]

    # Create containing folder if doesn't exist
    if not os.path.exists(orfs_fasta_path + "has_cog"):
        os.makedirs(orfs_fasta_path + "has_cog")
    
    # Save to file
    valid_orfs_path = orfs_fasta_path + "has_cog/" + sample + ".fna"
    SeqIO.write(valid_orfs, valid_orfs_path, "fasta")
    log("Saved ORFs with COGs for %s." % sample, t)

def extract_ORFs_without_COGs(sample):
    """ Extracts ORFs for genes that have COG assignment and saves their
    FASTA sequences. """
    sample = get_unique_sample(sample)
    sample_paths = get_sample_paths(sample)
    
    # Load sample data
    genes = load_sample_genes(sample)
    
    # Get ORFs with COGs
    t = time()
    has_cog = genes[genes.eggNOG_old == "unknown"].index.tolist()
    ORF_seqs = SeqIO.index(sample_paths["ORFs_fasta"], "fasta")
    valid_orfs = [ORF_seqs[orf] for orf in has_cog]

    # Create containing folder if doesn't exist
    if not os.path.exists(orfs_fasta_path + "no_cog_old"):
        os.makedirs(orfs_fasta_path + "no_cog_old")
    
    # Save to file
    valid_orfs_path = orfs_fasta_path + "no_cog_old/" + sample + ".fna"
    SeqIO.write(valid_orfs, valid_orfs_path, "fasta")
    log("Saved ORFs without COGs for %s." % sample, t)

def query_taxonomy(tax_ids, use_first_result=True):
    """ Returns the taxonomy for by querying the NCBI Taxonomy database. """
    
    # Make sure tax_ids is a list of strings
    if type(tax_ids) != list:
        tax_ids = [tax_ids]
    tax_ids = [str(int(tax_id)) for tax_id in tax_ids]
    
    # Query the taxonomy database via Entrez
    results = pd.DataFrame(Entrez.parse(Entrez.efetch("taxonomy", id=tax_ids)))
    
    # Resolve multiple results per query
    if use_first_result:
        results = results.groupby("TaxId").first().reset_index()
    
    # Return taxonomy tables
    taxonomy = [pd.DataFrame(lineage) for lineage in results.LineageEx]
    if len(taxonomy) == 1:
        taxonomy = taxonomy[0]
    return taxonomy

def update_genes_with_eggnog(sample, overwrite=False):
    """ Updates gene table for a given sample using BLAST output from querying
    the eggNOG database. """
        
    # Load current genes
    genes = load_sample_genes(sample)
    
    if "eggNOG_old" in genes and not overwrite:
        return
        
    # Load BLAST output
    blast_out = pd.read_table(blast_out_path + sample + ".tbl", names=blast_columns)
    
    # Get hits by lowest e-value
    keep_cols = ["hit", "e_value"]
    blast_hits = blast_out.groupby("query").first()[keep_cols]
    
    # Get COGs and functinal categories
    blast_hits = blast_hits.join(eggnog_members["#nog name"], on="hit")
    blast_hits = blast_hits.join(eggnog_funccat, on="#nog name")
    
    # Parse out eggNOG taxonomy ID
    blast_hits = blast_hits.join(blast_hits["hit"].str.extract("(?P<tax_id>\d+)\.(?P<gene_id>.*)").convert_objects(convert_numeric=True))
    
    # Get taxonomy for hits
    hits_tax = blast_hits.join(eggnog_tax, on="tax_id")
    hits_tax.rename(columns={"hit": "eggNOG_id", "e_value": "eggNOG_e_value", "#nog name": "eggNOG", "funccat": "eggNOG_funccat"}, inplace=True)
    hits_tax.drop(["tax_id", "gene_id"], axis=1, inplace=True)
    hits_tax.index.name = "gene_name"
    
    # Update with new taxonomy
    genes = genes.join(hits_tax, lsuffix="_old")
    
    # Save
    sample_genes_path = get_sample_paths(sample).genes
    genes.to_csv(sample_genes_path)
    
#%% Batch process samples
#samples = get_all_samples()
#samples = ["158256496-stool1", "158337416-stool1", "158337416-stool2", "158458797-stool1", "158479027-stool1", "158499257-stool1", "158499257-stool2", "158742018-stool1", "158802708-stool1", "158802708-stool2", "158883629-stool1", "158924089-stool1", "158944319-stool1", "158944319-stool2", "159005010-stool1", "159146620-stool1", "159166850-stool1", "159207311-stool1", "159207311-stool2", "159227541-stool1", "159227541-stool2", "159227541-stool3", "159247771-stool1", "159247771-stool2", "159268001-stool1", "159268001-stool2", "159369152-stool1", "159369152-stool2", "159490532-stool1", "159490532-stool2", "159510762-stool1", "159510762-stool2", "159551223-stool1", "159551223-stool2", "159571453-stool1", "159571453-stool2", "159591683-stool1", "159591683-stool2", "159611913-stool1", "159611913-stool2", "159632143-stool1", "159713063-stool1", "159733294-stool1", "159753524-stool1", "159753524-stool2", "159814214-stool1", "159814214-stool2", "159915365-stool1", "160158126-stool1", "160158126-stool2", "160178356-stool1", "160218816-stool1", "160319967-stool1", "160400887-stool1", "160421117-stool1", "160502038-stool1", "160582958-stool1", "160603188-stool1", "160643649-stool1", "160704339-stool1", "160765029-stool1", "246515023-stool1", "246515023-stool2", "338793263-stool1", "404239096-stool1", "508703490-stool1", "550534656-stool1", "604812005-stool1", "604812005-stool2", "638754422-stool1", "638754422-stool2", "675950834-stool1", "686765762-stool1", "686765762-stool2", "706846339-stool1", "737052003-stool1", "763435843-stool1", "763496533-stool1", "763496533-stool2", "763536994-stool1", "763536994-stool2", "763536994-stool3", "763577454-stool1", "763577454-stool2", "763597684-stool1", "763678604-stool1", "763678604-stool2", "763759525-stool1", "763759525-stool2", "763820215-stool1", "763820215-stool2", "763840445-stool1", "763840445-stool2", "763860675-stool1", "763860675-stool2", "763901136-stool1", "763901136-stool2", "763961826-stool1", "763961826-stool2", "763982056-stool1", "763982056-stool2", "764002286-stool1", "764042746-stool1", "764042746-stool2", "764062976-stool1", "764062976-stool2", "764143897-stool1", "764143897-stool2", "764184357-stool1", "764224817-stool1", "764224817-stool2", "764285508-stool1", "764325968-stool1", "764325968-stool2", "764447348-stool1", "764447348-stool2", "764487809-stool1", "764487809-stool2", "764508039-stool1", "764588959-stool1", "764649650-stool1", "764669880-stool1", "764669880-stool2", "764811490-stool1", "764892411-stool1", "765013792-stool1", "765034022-stool1", "765074482-stool1", "765074482-stool2", "765094712-stool1", "765094712-stool2", "765135172-stool1", "765560005-stool1", "765620695-stool1", "765640925-stool1", "765701615-stool1", "809635352-stool1", "823052294-stool1", "861967750-stool1"]
samples = [sample for sample in get_valid_samples() if "eggNOG_old" not in load_sample_genes(sample)]

failed = []; done = []
for sample in tqdm(samples):
    try:
        extract_ORFs_without_COGs(sample)
        done.append(sample)
    except:
        failed.append(sample)
        continue
    