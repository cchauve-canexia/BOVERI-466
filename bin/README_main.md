# The main.py script implements three tasks

1. latest_runs: computes the latest run directory for the given input runs

python path_to_main.py path_to_parameters_file path_to_log_file 'latest_runs'

It generates in the output directory specified in the parameters file, for each
MiSeq and NextsSeq input file F_manifest.csv (one per manifest, specified in the
parameters file), a file F_manifest_latest_run.csv

2. dump_files: generates a dump file of variants for each run in a pair of
matching MiSeq/NextSeq runs

python path_to_main.py path_to_parameters_file path_to_log_file 'dump_files'

Creates in the directory specified by DUMP_OUT in the parameters file a tsv dump
file for all variants of a run.
Output: csv file with 'run' column taht lists the run name and 'dump_file'
column that lists the path to the dump fle for the run.

3. variants: generates a csv file of all indel calls per manifest

python path_to_main.py path_to_parameters_file path_to_log_file 'variants'

Output is in the output directory specified in the parameters file.
All indel calls in a run belonging to a pair of matching runs are reported.

Output format: csv file where ach row is a variant seen in a sample f a pair of
matching sample, with columns
index: chromosome_position_reference_alternate
NextSeq_run: NextSeq run ID for the pair of matching runs where a sample contains the variant
MiSeq run: matching MiSeq run ID
NextSeq_sample: name of sample in the NextSeq run
MiSeq_sample: name of sample in the MiSeq run
  One of the two sample, maybe both, contain the variant
sample: sample ID independent of the considered run
gene: gene containing the variant
codon: codon affected by the variant
cDNA_change: cDNA change caused by the variant
vaf_n, score_n, coverage_n: VAF, Strelka score, coverage in the NextSeq sample (NaN if not called)
vaf_m, score_m, coverage_m: VAF, Strelka score, coverage in the MiSeq sample (NaN if not called)
