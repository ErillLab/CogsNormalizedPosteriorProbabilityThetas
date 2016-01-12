# CogsNormalizedPosteriorProbabilityThetas
Scripts for determining the normalized posterior probability of COGs above various thresholds.

## Initial steps
Make sure the samples have been scored with the IGC pipeline. The paths may need modifying.

Note: due to a rescoring and disk space issues, the data was split across 2 drives, and two paths were used.

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

