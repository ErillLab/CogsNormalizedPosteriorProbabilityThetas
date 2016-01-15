# CogsNormalizedPosteriorProbabilityThetas
Scripts for determining the normalized posterior probability of COGs above various thresholds.

## Initial steps
Make sure the samples have been scored with the IGC pipeline. The paths may need modifying.

Note: due to re-scoring and disk space issues, the data was split across 2 drives, and two paths were used.

The IGC pipeline needs to be run under Anaconda. There are convenience functions for running only "MH" samples or all the others. The pipeline can be run from an IDE like Spyder or simple python scripts that can be executed from the command line.

See also information about EggNOG Taxonomy below.

## Generate the data for the COG normalized posterior probabilities

Modify cogs_norm_post_prob_thetas.py to
1. Change the output directory as needed.
2. Set firmi true/false if Firmicutes or Gammaproteobacteria are desired.
3. Set skip to true/false if short promoters should/shouldn't be skipped.
4. Change what COGs will have their gene info saved ("out_cogs").
5. Change what COGs will have their data output to stdout for testing ("test_cogs").
6. The code currently processes MH-named samples. Remove the "not" in if-test ("if not sample[0:2]=="MH":") to skip these. Remove the if-test to process all samples.

The code is run by passing in the TF. Several predefined TFs have been set up in this script. New TFs should be set up in the same way.

Run under Anaconda.


## Summarize the COG normalized posterior probabilities

Once cogs_norm_post_prob_thetas.py has been run, the files in its output directory are processed with gather_norm_ratios_multis.py.

Modify this script to set mh_only to True or False as desired.

The script takes as input: the output directory (is the input for this script), a "file name root" to use to make the output of this script unique (becomes part of the filename),
the TF, and a COG to output for testing. Note that all COGs present input files are processed.

Run under Anaconda.


## Determine COG's genes and promoters

Edit the script run_get_prom_seqs.py to

1. Get the input directory (was the output directory for cogs_norm_post_prob_thetas.py)
2. Set the output directory for this script
3. Set the TF
4. Set the theta multiplier to grab (-6.0, -4.0, or -2.0 are good)
5. Set whether "firmi" or "gamma"

The script takes a COG or eggNOG as an argument.


## EggNOG Taxonomy

The script eggnog_taxonomy.py was used to integrate data obtained by BLASTing the eggNOG database to get COG/NOGs and taxonomy.

For the BLASTing to get the eggNOG hits we used DIAMOND (Buchfink, Xie & Huson, Nature Methods, 2015). Taxonomy was annotated in eggNOG v4 (Powell et al., 2013
<http://nar.oxfordjournals.org/content/42/D1/D231>).

The procedure was:

1. Build a DIAMOND database based on eggNOG v4.0 <http://eggnog.embl.de/version_4.0.beta/downloads.v4.html> (eggnogv4.proteins.all.fa.gz) which is a FASTA file of of AA seqs labelled with a tax_id.protein_id.

2. Extract ORFs with COGs from each sample and save them as separate (NT) FASTA files for processing.

3. Run DIAMOND (under a Ubuntu VM) using shell scripts on each of the samples.

4. Parse the output from DIAMOND (see eggnog_taxonomy.py), which were just the typical BLAST table output format, choose the best hit for each ORF in each sample and update the corresponding COG and taxonomy assignment in the genes table for each sample that was initially generated from the IGC data. The taxonomy assignment also required querying Entrez to get the different taxonomic levels from the tax_id that eggNOG provides (also done in eggnog_taxonomy.py).

