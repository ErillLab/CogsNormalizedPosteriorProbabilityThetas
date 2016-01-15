# -*- coding: utf-8 -*-
"""
This script automates processing of the IGC dataset.
Created on Tue Mar 17 19:13:08 2015
@author: Talmo
"""

import os
import sys
import time
import ast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from Bio import SeqIO
from PSSMScorer import PSSMScorer

#%% Configuration
# Output verbosity
verbosity = 1 # 0 = no output, 1 = minimal, 2 = debugging

# Base paths
base_path = "/home/cuda/2TB/metagenomics/" # metagenomics folder
#script_dir_path = os.path.dirname(os.path.abspath(__file__)) + os.path.sep
script_dir_path = base_path+"/mg_pipeline/" + os.path.sep
IGC_path = base_path + "IGC/"
alt_IGC_path = "/media/sf_D_DRIVE/metagenomics/IGC/" # 4TB drive

# Pipeline data paths
samples_index_path = script_dir_path + "data/samples_index.csv"
eggnog_tax_path = script_dir_path + "data/eggnogv4_taxonomy.csv"

# IGC Data paths
gene_summary_path = IGC_path + "3.IGC.AnnotationInfo/IGC.annotation.summary.v2"
scaftigs_path = IGC_path + "4.IndividualAssmeblies/"
orfs_fasta_path = IGC_path + "5.IndividualORFs/"

# Processed data paths
genes_path = IGC_path + "Genes/"
orfs_path = IGC_path + "ORFs/"
#operons_path = IGC_path + "Operons/"
operons_path = alt_IGC_path+"Operons/" ####################################################################

# PSSM scoring
binding_sites_path = base_path + "binding_sites/"
Firmicutes_LexA = PSSMScorer(binding_sites_path + "Firmicutes_LexA.txt")
Firmicutes_LexA.initialize_estimator(bg_mu=-17.681878, bg_sigma=8.2267625)
GammaProteobacteria_LexA = PSSMScorer(binding_sites_path + "GammaProteobacteria_LexA.txt")
GammaProteobacteria_LexA.initialize_estimator(bg_mu=-21.660313, bg_sigma=8.476820)
Grampos_LexA = PSSMScorer(binding_sites_path + "Grampositive_LexA.txt")
Firmicutes_Actino_LexA = PSSMScorer(binding_sites_path + "Firmicutes_Actinobacteria_LexA.txt")

YrkD_PSSM = PSSMScorer(binding_sites_path + "Bacillales_YrkD.txt", "YrkD")
YrkD_PSSM.initialize_estimator(bg_mu=-16.998722, bg_sigma=7.202773)

CsoR_PSSM = PSSMScorer(binding_sites_path + "Bacillales_CsoR.txt", "CsoR")
CsoR_PSSM.initialize_estimator(bg_mu=-16.362418, bg_sigma=6.652658)

CopY_PSSM = PSSMScorer(binding_sites_path + "CopY.txt", "CopY")
CopY_PSSM.initialize_estimator(bg_mu=-27.392480, bg_sigma=10.159172)

# Score results
scores_path = alt_IGC_path + "Scores/" 

# Operon prediction
threshold_IGI = 50 # max intergenic interval (bp)
promoter_region = (-250, +50) # bp relative to gene start

# BLASTing
blast_columns = ["query", "hit", "percentage_identity", "alignment_length", "num_mismatches", "num_gap_openings", "q.start", "q.end", "s.start", "s.end", "e_value", "score_bits"]

#%% General functions
# From the timeit module:
if sys.platform == 'win32':
    time = time.clock # On Windows, the best timer is time.clock
else:
    time = time.time # On most other platforms the best timer is time.time
    
def log(msg, start_time=None, verbosity_level=1):
    """ Convenience function for printing pipeline messages.
    
    Args:
        msg:
            the message string
        start_time:
            if provided, displays time elapsed since this time
        verbosity_level:
            if the verbosity is less than this, nothing will be printed
        
    Example:
        >>> start_time = time()
        >>> log("Done processing.", start_time)
        Done processing. [10.94s] """
        
    if verbosity_level <= verbosity:
        if start_time is not None:
            msg = msg + " [%.2fs]" % (time() - start_time)
        print msg

#%% Samples and paths
# Load pre-processed samples index
samples_index = pd.read_csv(samples_index_path, index_col="sample").dropna()
all_samples = samples_index.index.tolist()

# Build dictionary mapping all sample names and aliases to a unique sample name
sample_aliases = samples_index.index.to_series().to_dict()
sample_aliases.update(samples_index.alias.reset_index().set_index("alias").sample.to_dict())

def get_unique_sample_name(query_sample, error_not_found=False):
    """ Returns a unique sample name by resolving ambiguity with aliases.
    
    Args:
        query_sample: str
            Sample name or alias
        error_not_found: bool, default False
            If True, raises LookupError if query_sample is not found.
            
    Returns:
        unique_sample: str
            Unique and unambiguous sample name or None.
            This can be used to index into the samples_index table.
            
    Example:
        >>> get_unique_sample_name("MH192")
        'MH0192'
        >>> samples_index.loc[get_unique_sample_name("MH192")]
        ...
    """
    if query_sample in sample_aliases.keys():
        return sample_aliases[query_sample]
    elif error_not_found:
        raise LookupError("Sample not found: %s" % query_sample)
    else:
        return None

def get_sample_info(sample):
    """ Returns information about a sample.
    
    Args:
        sample: str
            Sample name or alias.
            
    Returns:
        sample_info: Series
            Sample information from samples_index table, e.g.:
                alias, SRA ID, BioSample ID, study, NCBI taxon ID/name, ...
    """
    sample_name = get_unique_sample_name(sample)
    sample_info = pd.Series({"sample": sample_name}).append(samples_index.loc[sample_name])
    return sample_info

def get_sample_paths(sample, original_paths=False):
    """ Returns paths to sample data files. 
    
    Args:
        sample: str
            Sample or alias. Raises an exception if not found.
        use_original: bool, default False
            If True, uses the original filenames from the downloaded files.
            If False, assumes ``standardize_paths()`` has renamed the files.
            
    Returns:
        paths: Series
            Absolute paths to the files associated with the sample. Keys:
            
            :scaftigs: Assembled scaftigs (in ``scaftigs_path``)
            :ORFs_fasta: ORF FASTA files (in ``orfs_fasta_path``)
            :genes: Summary of genes in sample (see ``save_sample_genes()``)
            :ORFs: ORF table without sequences (see ``parse_ORFs()``)
            :operons: Predicted operons (see ``predict_operons()``)
            
    Example:
        >>> get_sample_paths('MH192')
        scaftigs      /home/cuda/2TB/metagenomics/IGC/4.IndividualAs...
        ORFs_fasta    /home/cuda/2TB/metagenomics/IGC/5.IndividualOR...
        genes          /home/cuda/2TB/metagenomics/IGC/Genes/MH0192.csv
        ORFs            /home/cuda/2TB/metagenomics/IGC/ORFs/MH0192.csv
        operons       /home/cuda/2TB/metagenomics/IGC/Operons/MH0192...
        Name: MH0192, dtype: object
    """
    # Get unique sample name
    sample = get_unique_sample_name(sample, True)
    
    # Initialize paths Series
    paths = pd.Series(name=sample)
    
    # IGC data
    if original_paths:
        paths["scaftigs"] = scaftigs_path + samples_index.at[sample, "Assemblies_filenames"]
        paths["ORFs_fasta"] = orfs_fasta_path + samples_index.at[sample, "ORFs_filenames"]
    else:
        paths["scaftigs"] = scaftigs_path + sample + ".fna"
        paths["ORFs_fasta"] = orfs_fasta_path + sample + ".fna"
    
    # Processed data
    paths["genes"] = genes_path + sample + ".csv"
    paths["ORFs"] = orfs_path + sample + ".csv"
    paths["operons"] = operons_path + sample + ".csv"
    
    # Make sure paths are absolute
    paths = paths.apply(os.path.abspath)
    
    return paths

def standardize_paths(sample):
    """ Renames the data files for a sample to their unique sample name.
    
    Affects scaftigs (assemblies) and ORF sequence FASTA files.
    
    Args:
        sample: str
            Sample or alias
    """
    # Get paths
    original_paths = get_sample_paths(sample, original_paths=True)
    new_paths = get_sample_paths(sample, original_paths=False)
    
    # Rename scaftigs
    if not os.path.exists(new_paths["scaftigs"]):
        if os.path.exists(original_paths["scaftigs"]):
            try:
                os.rename(original_paths["scaftigs"], new_paths["scaftigs"])
                log("Renamed assembled scaftig for sample %s." % sample, verbosity_level=2)
            except:
                log("Could not rename original assembled scaftig for sample %s." % sample, verbosity_level=1)
        else:
            log("Could not find original assembled scaftig for sample %s." % sample, verbosity_level=1)
    
    # Rename ORF sequences
    if not os.path.exists(new_paths["ORFs_fasta"]):
        if os.path.exists(original_paths["ORFs_fasta"]):
            try:
                os.rename(original_paths["ORFs_fasta"], new_paths["ORFs_fasta"])
                log("Renamed ORF sequences for sample %s." % sample, verbosity_level=2)
            except:
                log("Could not rename original ORF sequences for sample %s." % sample, verbosity_level=1)
        else:
            log("Could not find original ORF sequences for sample %s." % sample, verbosity_level=1)


def has_data(sample):
    """ Checks if the sample has its data files. 
    
    Args:
        sample: str
            The sample to be queried
            
    Returns:
        data_exists: Series
            Same keys as returned by ``get_sample_paths()`` but will be True
            for the ones for which the file exists.
    """
    return get_sample_paths(sample).apply(os.path.exists)

def check_samples_data():
    """ Returns a table of all samples and which data files are present. """
    return pd.Series(all_samples, index=all_samples).apply(has_data)

#def get_valid_samples():
#    """ Returns a list of all samples that have all data files present. """
#    valid = check_samples_data().all(1)
#    return valid.index[valid].tolist()

def get_MetaHit(study=2010):
    """ Returns the names of the samples from the MetaHit database. 
    
    Args:
        study: int or str, default 2010
            Year of the study to return. Studies:
            
            :2010: Qin et al. (2010) | doi:10.1038/nature08821
            :2012: Qin et al. (2012) | doi:10.1038/nature11450
            :2013: Le Chatelier et al. (2013) | doi:10.1038/nature12506
            :"all": Returns all studies.
    
    Returns:
        MetaHit_samples: list
            The samples corresponding to the study selected
    """
    
    queries = {"2010": "A human gut microbial gene catalogue", \
               "2012": "A metagenome-wide association", \
               "2013": "Richness of human gut"}
    
    if study == 'all':
        return [item for inner_list in [get_MetaHit(k) for k in queries.keys()] for item in inner_list]
    
    return samples_index[samples_index.study.str.contains(queries[str(study)])].index.tolist()

def get_all_samples(HMP=True):
    """ Returns a list of all samples with data. 
    
    Args:
        HMP: bool, default True
            If False, filters out the samples that belong to the Human
            Microbiome Project.
    """
    if not HMP:
        from_HMP = samples_index.study.str.contains("Human Microbiome Project")
        return samples_index.index[~from_HMP].index.tolist()
    return all_samples

#def standardize_paths(sample):

#%% Genes processing
# Regular expression that matches all sample names and aliases
samples_regex = "(" + "|".join(sample_aliases.keys()) + ")"
def load_gene_summary(limit=None, extract_samples=False):
    """ Loads the integrated gene summary.
    
    Warning: This may take a while as there are 9879896 rows in the table.
    
    Args:
        limit: int, default None
            The number of rows to return. If None, returns all.
        extract_sample: bool, default False
            Tries to extract the sample that each gene comes from and appends 
            it to the "sample" column (see ``extract_sample()``).
    
    Returns:
        gene_summary: DataFrame
            The IGC gene summaries table indexed by gene_id.
        
    Example:
        >>> print load_gene_summary(1)
    ================  ===========  ============  =============  =======  =======  ======  =======  ==============  ===============  ================  =======================================================================  ================
    gene_name         gene_length  completeness  cohort_origin  phylum   genus    kegg    eggNOG   sample_freq     individual_freq  eggNOG_funccat    kegg_funccat                                                             cohort_assembled
    ================  ===========  ============  =============  =======  =======  ======  =======  ==============  ===============  ================  =======================================================================  ================
    T2D-6A_GL0083352  88230        Complete      CHN            unknown  unknown  K01824  COG5184  0.224151539069  0.236448598131   Lipid Metabolism  Cell cycle control, cell division, chromosome partitioning;Cytoskeleton  EUR;CHN;USA
    ================  ===========  ============  =============  =======  =======  ======  =======  ==============  ===============  ================  =======================================================================  ================
    """
    
    t = time()
    gene_summary_cols = ['gene_id', 'gene_name', 'gene_length', 'completeness', \
        'cohort_origin', 'phylum', 'genus', 'kegg', 'eggNOG', 'sample_freq', \
        'individual_freq', 'eggNOG_funccat', 'kegg_funccat', 'cohort_assembled']
    gene_summary = pd.read_table(gene_summary_path, names=gene_summary_cols, index_col='gene_id', nrows=limit)
    log("Loaded %d genes from summary." % (gene_summary.shape[0]), t)
    
    if extract_samples:
        gene_summary["sample"] = extract_sample(gene_summary["gene_name"])
    
    return gene_summary

def extract_sample(gene_names):
    """ Returns samples given one or more gene_names. 
    
    Args:
        gene_names: str, list, dict, Series
            One or more gene_names
            
    Returns:
        samples: str, list, dict, Series
            Returns a Series if input was not of above types.
            
    Example:
        >>> genes = load_sample_genes('MH192')
        >>> extract_sample(genes.index)
        0     MH0192
        ...
        21241    MH0192
        Name: gene_name, Length: 21242, dtype: object
        >>> extract_sample(['MH0192_GL0000004', 'MH0192_GL0000005', 'MH0193_GL0000001'])
        ['MH0192', 'MH0192', 'MH0193']
    """
    
    # Convert to Series
    if type(gene_names) != pd.Series:
        original_type = type(gene_names)
        gene_names = pd.Series(gene_names)
    
    # Match against regex
    t = time()
    samples = gene_names.str.extract(samples_regex)
    v = [1, 2][(time() - t) < 1.0]
    log("Extracted sample regex.", t, v)
    
    # Get unique names for samples
    t = time()
    samples = samples.apply(get_unique_sample_name)
    v = [1, 2][(time() - t) < 1.0]
    log("Disambiguated sample names.", t, v)
    
    if original_type == str:
        return samples.ix[0]
    elif original_type == list:
        return samples.tolist()
    elif original_type == dict:
        return samples.to_dict()
    else:
        return samples


def save_sample_genes(overwrite=False, gene_summary=None):
    """ Splits the integrated gene summary into individual tables for each 
    sample.
    
    To load these files, see ``load_sample_genes()``.
    
    See ``genes_path`` for path to data files.
    
    Args:
        overwrite: bool, default False
            If True, will overwrite pre-existing gene files.
        gene_summary: DataFrame, default None
            Accepts a pre-loaded gene_summary DataFrame (see 
            ``load_gene_summary()``). If None, calls that function to load it.
    
    """
    
    # Create parent folder if it does not exist
    if not os.path.exists(genes_path):
        os.makedirs(genes_path)
    
    # Load the combined gene summary table
    if gene_summary is None:
        gene_summary = load_gene_summary()

    # Extract samples from gene names
    if "sample" not in gene_summary.columns:
        gene_summary["sample"] = extract_sample(gene_summary["gene_name"])

    # Drop genes without sample
    n_total = len(gene_summary)
    gene_summary = gene_summary[~gene_summary.sample.isnull()]
    n_after = len(gene_summary)
    log("%d/%d genes without sample (%.2f%%)." % (n_total - n_after, n_total, 100.0*(n_total - n_after)/n_total), None, 2)
    
    empty_samples = samples_index[~(samples_index.Assemblies_filenames.isnull() | samples_index.ORFs_filenames.isnull())].index.diff(gene_summary.sample.unique())
    log("%d/%d samples without genes (%.2f%%)." % (len(empty_samples), len(samples_index), 100.0*len(empty_samples)/len(samples_index)), None, 2)
    
    # Group by sample and save data
    t = time()
    gene_summary = gene_summary.groupby("sample", sort=False)
    for sample, sample_genes in gene_summary:
        sample_genes_path = get_sample_paths(sample).genes
        
        if overwrite or not os.path.exists(sample_genes_path):
            sample_genes.to_csv(sample_genes_path)
            
    log("Saved genes for %d samples individually." % gene_summary.ngroups, t)

def load_sample_genes(sample):
    """ Loads the genes table for a sample.
    
    See ``save_sample_genes()`` for table creation.
    
    Args:
        sample: str
            Sample name or alias
        
    Returns:
        genes: DataFrame
            Listing and information about each called gene in the sample.
     """
    # Get path
    sample_genes_path = get_sample_paths(sample).genes
    
    # Check if it exists
    if not os.path.exists(sample_genes_path):
        raise EnvironmentError("Genes file for %s does not exist: %s" % (sample, sample_genes_path))
        
    # Load and return
    sample_genes = pd.read_csv(get_sample_paths(sample).genes, index_col="gene_name")
    sample_genes.sort_index(inplace=True)
    return sample_genes

def get_genes(sample):
    """ Loads the genes table for a sample. Shortcut for
    ``load_sample_genes()``.
    
    See ``save_sample_genes()`` for table creation.
    
    Args:
        sample: str
            Sample name or alias
        
    Returns:
        genes: DataFrame
            Listing and information about each called gene in the sample.
    """
    return load_sample_genes(sample)

#%% Scaftigs
def get_scaftigs(sample, index_only=False, convert_to_str=True):
    """ Loads the assembled scaftig nucleotide sequences.
    
    Args:
        sample: str
            Sample name or alias
        index_only: bool, default False
            If True, will read in the whole file and parse it into a dict.
            If False, will only index the file. This uses less memory by not
            loading sequences into memory (see Bio.SeqIO.index).
        convert_to_str: bool, default True
            Converts Bio.Seq objects to strings. This improves performance
            by avoiding the overhead of Bio.Seq objects at the cost of utility
            functions in Bio.Seq. Only applicable when index_only is False.
    
    Returns:
        scaftigs: dict or dict-like
            The scaftig sequences indexed by their name.
            Returns None if not found.
    """
    t = time()
    
    # Get unique sample and paths
    sample = get_unique_sample_name(sample)
    file_path = get_sample_paths(sample)["scaftigs"]
    
    # Check if file exists
    if not os.path.exists(file_path):
        log("Could not find scaftigs for sample %s in: %s" % (sample, file_path), verbosity_level=1)
        return
        
    if index_only:
        # Create file index
        scaftigs = SeqIO.index(file_path, "fasta")
        log("Indexed scaftigs for sample %s." % sample, t, 2)
    else:
        # Read and parse file
        scaftigs = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
        
        # Downcast to strings
        if convert_to_str:
            for scaf, seq in scaftigs.items():
                scaftigs[scaf] = str(seq.seq).upper()
        
        log("Read and parsed scaftigs for sample %s." % sample, t, 2)
    
    return scaftigs

#%% ORFs
def get_ORF_seqs(sample, index_only=False, convert_to_str=True):
    """ Loads the ORF FASTA nucleotide sequences.
    
    Args:
        sample: str
            Sample name or alias
        index_only: bool, default False
            If True, will read in the whole file and parse it into a dict.
            If False, will only index the file. This uses less memory by not
            loading sequences into memory (see Bio.SeqIO.index).
        convert_to_str: bool, default True
            Converts Bio.Seq objects to strings. This improves performance
            by avoiding the overhead of Bio.Seq objects at the cost of utility
            functions in Bio.Seq. Only applicable when index_only is False.
    
    Returns:
        ORF_seqs: dict or dict-like
            The ORF sequences indexed by their name in the ORFs table.
            Returns None if not found.
    """
    t = time()
    
    # Get unique sample and paths
    sample = get_unique_sample_name(sample)
    file_path = get_sample_paths(sample)["ORFs_fasta"]
    
    # Check if file exists
    if not os.path.exists(file_path):
        log("Could not find ORF sequences for sample %s in: %s" % (sample, file_path), verbosity_level=1)
        return
        
    if index_only:
        # Create file index
        ORF_seqs = SeqIO.index(file_path, "fasta")
        log("Indexed ORF sequences for sample %s." % sample, t, 2)
    else:
        # Read and parse file
        ORF_seqs = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
        
        # Downcast to strings
        if convert_to_str:
            for orf, seq in ORF_seqs.items():
                ORF_seqs[orf] = str(seq.seq).upper()
        
        log("Read and parsed ORF sequences for sample %s." % sample, t, 2)
    
    return ORF_seqs

def parse_ORFs(sample, overwrite=False):
    """ Parses the FASTA headers in the ORF FASTA file and saves the data.
    
    Args:
        sample: str
            Sample name or alias
        overwrite: boolean, default False
            If True, will overwrite saved ORF table, otherwise will attempt to
            load the table from cache.
            
    Returns:
        ORFs: DataFrame
            Table with data on the open reading frames and the columns:
                gene_name, orf_type, scaffold, start, end, completeness
                
    Example:
        >>> ORFs = parse_ORFs('MH0001')
        Parsed 21452 ORFs from MH0001. [0.44s]
        >>> ORFs.head(1)
        ...
        """
    t = time()
    
    # Get unique sample and paths
    sample = get_unique_sample_name(sample)
    sample_paths = get_sample_paths(sample)
    
    # Check for cache
    if not overwrite and os.path.exists(sample_paths["ORFs"]):
        ORFs = pd.read_csv(sample_paths["ORFs"], index_col=0)
        log("Loaded %d cached ORFs from %s." % (ORFs.shape[0], sample), t)
        return ORFs
    
    # Read in FASTA headers
    with open(sample_paths["ORFs_fasta"]) as f:
        lines = [line[1:-1] for line in f if line[0] == '>']
    
    # Check header format
    ORF_pattern = '(?P<gene_name>[^ ]+)[ ]+\\[(?P<orf_type>[^\]]+)\\][ ]+locus=(?P<scaffold>[^:]+):(?P<start>[^:]+):(?P<end>[^:]+):(?P<strand>[^:\[]+)\\[(?P<completeness>[^\[\]]+)'
    if "]" not in lines[0]:
        ORF_pattern = '(?P<gene_name>[^\t]+)\tstrand:(?P<strand>[+-]) start:(?P<start>\d+) stop:(?P<end>\d+) length:(?P<length>\d+) start_codon:(?P<start_codon>\w+) stop_codon:(?P<stop_codon>\w+) gene_type:(?P<gene_type>\w+)'
    
    # Parse FASTA headers for ORFs
    ORFs = pd.Series(lines).str.extract(ORF_pattern).convert_objects(convert_numeric=True)
    
    # Standardize alternative format
    if "]" not in lines[0]:
        # completeness
        ORFs.loc[(ORFs.start_codon == 'no') & (ORFs.stop_codon == 'no'), 'completeness'] = "Lack both ends"
        ORFs.loc[(ORFs.start_codon == 'yes') & (ORFs.stop_codon == 'no'), 'completeness'] = "Lack 3'-end"
        ORFs.loc[(ORFs.start_codon == 'no') & (ORFs.stop_codon == 'yes'), 'completeness'] = "Lack 5'-end"
        ORFs.loc[(ORFs.start_codon == 'yes') & (ORFs.stop_codon == 'yes'), 'completeness'] = "Complete"
        
        # orf_type
        ORFs.loc[:, "orf_type"] = "gene"
        
        # scaffold
        ORFs.loc[:, "scaffold"] = ORFs.gene_name.str.extract("(?P<gene_name>.+)_gene")
        
        # Re-order columns
        ORFs = ORFs.loc[:, ['gene_name', 'orf_type', 'scaffold', 'start', 'end', 'strand', 'completeness']]
    
    # Set gene_name as index
    ORFs.set_index("gene_name", inplace=True)
    
    # Create parent folder if it does not exist
    if not os.path.exists(orfs_path):
        os.makedirs(orfs_path)
    
    # Save
    ORFs.to_csv(sample_paths["ORFs"])
    log("Parsed %d ORFs from %s." % (ORFs.shape[0], sample), t)
    
    return ORFs
    
def get_ORFs(sample):
    """ Loads the ORFs table for a sample. See ``parse_orfs()`` to create an
    ORFs table for a sample.
    
    Args:
        sample: str
            Sample name or alias
            
    Returns:
        ORFs: DataFrame
            Table with a listing and information about each ORF predicted in
            the sample.
    """
    t = time()
    
    # Get unique sample and paths
    sample = get_unique_sample_name(sample)
    sample_paths = get_sample_paths(sample)
    
    # Check for cache
    if os.path.exists(sample_paths["ORFs"]):
        ORFs = pd.read_csv(sample_paths["ORFs"], index_col=0)
        log("Loaded %d cached ORFs from %s." % (ORFs.shape[0], sample), t, 2)
        return ORFs

#%% Operons
def predict_operons(sample, overwrite=False):
    """ Predicts operons for the given sample. 
    
    Args:
        sample: str
            Sample for which to predict operons.
        overwrite: bool, default False
            If False, loads predictions from disk if they exist.
    """
    t = time()
    
    # Get paths
    sample_paths = get_sample_paths(sample)
    
    # Check for cache
    if not overwrite and os.path.exists(sample_paths["operons"]):
        operons = pd.read_csv(sample_paths["operons"], index_col=0, converters={'genes':ast.literal_eval})
        log("Loaded %d cached operons from %s." % (operons.shape[0], sample), t)
        return operons
    
    # Load sample scaffolds (scafftigs)
    scaffolds = SeqIO.index(sample_paths["scaftigs"], 'fasta')
    
    # Load the predicted open reading frames
    ORFs = parse_ORFs(sample)
    
    # Group by 'scaffold' column and do groupwise operations
    ORFs = ORFs.groupby('scaffold', sort=False)
    n_ORFs = ORFs.ngroups
    print "Number of ORFs", n_ORFs
    ocnt = 0
    scnt = 0
    ofname = operons_path+ sample+".csv"
    print "Creating", ofname
    out_fp = open(ofname, "w")
    out_fp.write("operon,genes,strand,promoter_start,promoter_end,promoter_seq,head_completeness\n")
    #operon_tables = []; i = 0
    for scaffold_name, scaf_ORFs in ORFs:
        scnt += 1
        #print "Process scaffold", scnt, scaffold_name
        scaf_info = scaffolds[scaffold_name]
        scaf_len = len(scaf_info)
        for strand, strand_ORFs in scaf_ORFs.groupby('strand', sort=False):
            # Sort by start values
            if strand == "-":
                strand_ORFs.sort('end', ascending=False, inplace=True)
                # Compute intergenic intervals
                IGIs = strand_ORFs.start[:-1].values - strand_ORFs.end[1:].values 
            else:
                strand_ORFs.sort('start', inplace=True) 
                # Compute intergenic intervals
                IGIs = strand_ORFs.start[1:].values - strand_ORFs.end[:-1].values 
        
            # Find operon gene indices
            head_genes = np.where(IGIs > threshold_IGI)[0] + 1
            operons_start = np.hstack(([0], head_genes))
            operons_end = np.hstack((head_genes, [strand_ORFs.shape[0]])) #np.hstack((head_genes, [scaf_ORFs.shape[0]]))
            if strand == "-":
                head_gene_start = strand_ORFs.end.iloc[operons_start]
            else:
                head_gene_start = strand_ORFs.start.iloc[operons_start] #scaf_ORFs.start.iloc[operons_start]
            if scaffold_name == "scaffold42382_1":
              print scaffold_name, strand, head_gene_start
            # Get promoter regions
            rmult = 1
            if strand == "-":
                rmult = -1 # Go the other way for reverse strands
            promoter_start = (head_gene_start + rmult*promoter_region[0]).clip(1).reset_index(drop=True)
            promoter_end = (head_gene_start + rmult*promoter_region[1]).reset_index(drop=True)
            if strand == "-":
                for i in range(0, len(promoter_start)):
                    if promoter_start[i] > scaf_len:
                        promoter_start[i] = scaf_len
                    # Talmo indicated that if the operon has start>end, to swap them. So 1795,1710 becomes 1710,1795
                    #promoter_seq = [str(scaffolds[scaffold_name][s-1:e].seq) for s, e in zip(promoter_end, promoter_start)]
                    temp = promoter_start[i]
                    promoter_start[i] = promoter_end[i]
                    promoter_end[i] = temp
            else:
                for i in range(0, len(promoter_end)):
                    if promoter_end[i] > scaf_len:
                        promoter_end[i] = scaf_len
            promoter_seq = [str(scaffolds[scaffold_name][s-1:e].seq) for s, e in zip(promoter_start, promoter_end)]
        
            # Build operon table
            strand_operons = pd.DataFrame()
            strand_operons['genes'] = [strand_ORFs.index[op_start:op_end].tolist() for op_start, op_end in zip(operons_start, operons_end)]
            strand_operons['strand'] = strand
            strand_operons['promoter_start'] = promoter_start
            strand_operons['promoter_end'] = promoter_end
            strand_operons['promoter_seq'] = promoter_seq
            strand_operons['head_completeness'] = strand_ORFs.ix[op_start, 'completeness']
            

            slen = len(strand_operons)
            for i in range(0, slen):
                glen = len(strand_operons['genes'][i])
                gout_str = "["
                for g in range(0,glen):
                    gout_str += "'"+strand_operons['genes'][i][g]+"', "
                gout_str = gout_str[:-2] + "]"
                if glen > 1:
                    gout_str = '"'+gout_str+'"'
                out_str=str(ocnt)+","+gout_str+","+strand+","+str(promoter_start[i])+","+str(promoter_end[i])+","+promoter_seq[i]+","+strand_operons['head_completeness'][i]
                out_fp.write("%s\n" % out_str)
                ocnt+=1
            # Append to list of operon tables
            #operon_tables.append(strand_operons)
        
        #i = i + 1
        #if verbosity >= 2 and i % (n_ORFs / 5) == 0:
        #    log("%.2f%%: %d/%d" % (float(i) / n_ORFs * 100, i, n_ORFs), None, 2)
        
    # Merge operon tables
    #operons = pd.concat(operon_tables, ignore_index=True)
    #operons.index.name = "operon"
    
    # Create parent folder if it does not exist
    #if not os.path.exists(operons_path):
    #    os.makedirs(operons_path)
        
    # Save
    #operons.to_csv(sample_paths["operons"])
    
    #log("Predicted %d operons for %s." % (operons.shape[0], sample), t)
    #return operons
    out_fp.close()
    return []
    
def get_operons(sample):
    """ Returns the operons for a sample if it exists. 
    
    Args:
        sample: str
            Sample name or alias
            
    Returns:
        operons: DataFrame
            Table with listing of predicted operons for a sample.
    """
    
    # Get path
    operons_path = get_sample_paths(sample)["operons"]
    
    # Check for cache
    if os.path.exists(operons_path):
        t = time()
        operons = pd.read_csv(operons_path, index_col=0, converters={'genes':ast.literal_eval})
        log("Loaded %d cached operons from %s." % (operons.shape[0], sample), t, 2)
        return operons
    else:
        return predict_operons(sample)

def get_genes2operon(genes, operons):
    """ Returns a Series that maps gene names to the operons they belong to. 
    
    Args:
        genes: DataFrame
            See ``get_genes()``.
        operons: DataFrame
            See ``get_operons()``.
        
    Returns:
        genes2operon: Series
            Operon numbers indexed by the gene names.
        """
    return pd.Series(index=np.hstack(operons.genes), data=np.repeat(operons.index.values, operons.genes.map(len)))

#%% BLASTing
def get_bash_array(samples):
    """ Returns a bash array of the samples.
    
    Args:
        samples: list
            List of sample names.
    
    Returns:
        bash_array: str
            A formatted bash array containing all samples names as strings."""
    return 'samples=("' + '" "'.join(samples) + '")'

#%% PSSM scoring
def score_sample(sample, PSSM, overwrite=False, soft_max=False):
    """ Scores the promoters of all the operons in a sample using the PSSM. 
    
    Args:
        sample: str
            Name of the sample to score
        PSSM: PSSMScorer
            PSSM to use for scoring
        overwrite: bool, default False
            If False, will load cached scores from disk
            
    Returns:
        scores: DataFrame
            The scores for the sample
        """
        
    # Validate sample name
    sample = get_unique_sample_name(sample)
    
    # Get data file path
    pssm_scores_path = scores_path + PSSM.name + "/"
    sample_scores_path = pssm_scores_path + sample + ".h5"
    
    # Check for cache
    if not overwrite and os.path.exists(sample_scores_path):
        print "Have scores"
        return get_sample_scores(sample, PSSM)
        
    # Load operons for sample
    operons = get_operons(sample)
    
    # Extract promoters and capitalize
    seqs = operons.promoter_seq.str.upper()
    
    # Score all the promoter sequences
    t = time()
    scores = seqs.apply(PSSM.score).apply(pd.Series, args=([["+","-"]]))
    log("Scored %s: %d sequences, %d bp." % (sample, len(seqs), seqs.apply(len).sum()), t, 1)
    
    # Compute soft-max
    if soft_max:
        scores["soft_max"] = scores.applymap(np.exp).sum(1).map(np.log)
    
    # Check if containing folder exists
    if not os.path.exists(pssm_scores_path):
        os.makedirs(pssm_scores_path)
    
    # Check if scores file exists
    if os.path.exists(sample_scores_path):
        os.remove(sample_scores_path) # delete existing
        
    # Save to HDF5
    t = time()
    scores.to_hdf(sample_scores_path, "table", format="fixed", append=False, mode="w")
    log("Saved scores for %s: %s" % (sample, sample_scores_path), t, 1)
    
    return scores
    
def get_scores_path(sample, PSSM_name):
    """ Returns the path to the scores for the given sample and PSSM name. """
    
    # Get name from PSSM instance
    if isinstance(PSSM_name, PSSMScorer):
        PSSM_name = PSSM_name.name
    
    # Validate sample name
    sample = get_unique_sample_name(sample)
    
    # Get data file path
    pssm_scores_path = scores_path + PSSM_name + "/"
    sample_scores_path = pssm_scores_path + sample + ".h5"
    
    return sample_scores_path
    
def get_sample_scores(sample, PSSM_name, soft_max=True):
    """ Loads scores for a given sample and PSSM name. """
    
    # Get path to scores
    sample_scores_path = get_scores_path(sample, PSSM_name)
    
    # Check if file exists
    if not has_score(sample, PSSM_name):
        raise EnvironmentError("Scores could not be found.")
        
    # Load scores
    t = time()
    scores = pd.read_hdf(sample_scores_path, "table")
    
    # Take soft max
    if soft_max:
        scores = scores.applymap(np.exp).sum(1).map(np.log)
    
    log("Loaded cached %s scores for %s." % (PSSM_name, sample), t, 2)
    
    return scores
    
def has_score(sample, PSSM_name):
    """ Returns True if the sample has scores under a PSSM name. """
    # Check if file exists
    return os.path.exists(get_scores_path(sample, PSSM_name))

def get_all_with_scores(PSSM_name):
    """ Returns a list with the names of all samples that have a score for the
    given PSSM name. """
    return [sample for sample in get_all_samples() if has_score(sample, PSSM_name)]

#%% Summary stats
def get_sample_summary(sample, PSSM):
    sample = get_unique_sample_name(sample)
    
    # Load sample data
    genes = load_sample_genes(sample)
    ORFs = parse_ORFs(sample)
    operons = predict_operons(sample)
    scores = get_sample_scores(sample, PSSM)
    
    stats = pd.Series()
    stats["operons"] = len(operons)
    stats["ORFs"] = len(ORFs)
    stats["genes"] = len(genes)
    
    stats["has_phylum"] = (genes.phylum != "unknown").sum()
    stats["has_COG"] = (genes.eggNOG != "unknown").sum()
    stats["has_both"] = ((genes.phylum != "unknown") & (genes.eggNOG != "unknown")).sum()
    
    # Find hits
    hit_threshold = 8.0
    hits = scores.applymap(lambda x: x >= hit_threshold)
    stats["hit_sites"] = hits.applymap(sum).sum(1).sum()
    
    hit_operons = hits.applymap(sum).sum(1) > 0
    stats["hit_operons"] = hit_operons.sum()
    
    hit_ORFs = np.hstack(operons.loc[hit_operons, "genes"])
    stats["hit_ORFs"] = len(hit_ORFs)
    
    hit_genes = genes.index.intersection(hit_ORFs)
    stats["hit_genes"] = len(hit_genes)
    
    stats["hit_has_phylum"] = (genes.loc[hit_genes, "phylum"] != "unknown").sum()
    stats["hit_has_COG"] = (genes.loc[hit_genes, "eggNOG"] != "unknown").sum()
    stats["hit_has_both"] = ((genes.loc[hit_genes, "eggNOG"] != "unknown") & (genes.loc[hit_genes, "phylum"] != "unknown")).sum()
    
    return stats

def get_samples_summary(samples="all", PSSM=Firmicutes_LexA):
    if samples == "all":
        samples = get_all_samples()
        
    stats = pd.Series(samples).apply(lambda x: get_sample_summary(x, PSSM))
    stats.index = samples
    stats.index.name = "sample"
    
    return stats
 
def do_score_test():
    sample_name="MH0239"
    bsite_pssm = GammaProteobacteria_LexA # Check scores_path near top
    bsite_pssm = Firmicutes_LexA
    rescore = True
    scores = score_sample(sample_name, bsite_pssm, overwrite=rescore)
    #sm_scores = get_sample_scores(sample_name, bsite, soft_max=True)

def run_predict_operons(start_at, stop_before):
  samples = get_all_samples()
  for sample in samples:
    if not sample[0:2]=="MH":
      continue
    snum = int(sample[2:])
    if (snum>=start_at and snum<stop_before):
      print sample
      the_ops = predict_operons(sample, False)
def run_nonMH_predict_operons():
  samples = get_all_samples()
  for sample in samples:
    if sample[0:2]=="MH":
      continue
    print sample
    the_ops = predict_operons(sample, False)
      
def run_scores(start_at, stop_before, TF, bsite_fname):
  print "binding sites path:", binding_sites_path, "fname:", bsite_fname
  bsite_pssm = PSSMScorer(binding_sites_path + bsite_fname, TF)
  samples = get_all_samples()
  for sample in samples:
    if  not sample[0:2]=="MH": 
      continue
    snum = int(sample[2:])
    if (snum>=start_at and snum<stop_before):
      print sample
      the_scores = score_sample(sample, bsite_pssm, False)
      
def run_nonMH_scores(TF, bsite_fname):
  print "binding sites path:", binding_sites_path, "fname:", bsite_fname
  bsite_pssm = PSSMScorer(binding_sites_path + bsite_fname, TF)
  samples = get_all_samples()
  for sample in samples:
    if sample[0:2]=="MH": # 
      continue
    print sample
    the_scores = score_sample(sample, bsite_pssm, False)
      
      
def change_score_path(which_drive):
    global scores_path
    if which_drive == "4TB":
        scores_path = "/media/sf_D_DRIVE/metagenomics/IGC/Scores/"
    else:
        scores_path = IGC_path + "Scores/"
