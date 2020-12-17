'''
Collecting data to generate the indel results of strelka on a list of  MiSeq runs.
Identifying the dump variants file and, if none, the output and raw data files
for each  runs.
Computing missing dump files.
Generating a CSV file of unfiltered indel calls
'''

# Standard imports
from collections import defaultdict
import os
import subprocess
import sys
import tempfile

# Third-party imports
import pandas as pd
import numpy as np

# Local imports
import utils
from data_format import *

# Computing the latest runs files for all runs of matching pairs
def compute_latest_runs(parameters, log):
    MANIFESTS = list(utils.get_manifests(parameters).keys())
    # Data files, indexed by manifest
    MISEQ_RUNS = utils.get_MiSeq_runs_file(parameters)
    MISEQ_EXCLUDED_RUNS = utils.get_MiSeq_excluded_runs_file(parameters)
    # Results files, indexed by manifest
    MISEQ_LATEST_RUNS = utils.get_MiSeq_latest_runs_file(parameters)
    for manifest in MANIFESTS:
        # Create a temporary non-excluded run file
        NON_EXCLUDED_RUNS = tempfile.NamedTemporaryFile(mode='w', delete=False)
        EXCLUDED_RUNS = []
        for _, run_data in pd.read_csv(MISEQ_EXCLUDED_RUNS[manifest], sep=',').iterrows():
            EXCLUDED_RUNS.append(run_data[EXC_RUN])
            runs_to_remove = []
        MISEQ_RUNS_DF =  pd.read_csv(MISEQ_RUNS[manifest], sep=',')
        for idx, run_data in MISEQ_RUNS_DF.iterrows():
            if run_data[RUN_ID]  in EXCLUDED_RUNS:
                runs_to_remove.append(idx)
        FILTERED_MISEQ_RUNS_DF = MISEQ_RUNS_DF.drop(runs_to_remove, axis=0)
        FILTERED_MISEQ_RUNS_DF.to_csv(NON_EXCLUDED_RUNS, index=False, columns=FILTERED_MISEQ_RUNS_DF.columns)
        NON_EXCLUDED_RUNS.close()
        # MiSeq run
        miseq_in, miseq_out = NON_EXCLUDED_RUNS.name, MISEQ_LATEST_RUNS[manifest]
        latest_runs_cmd = utils.get_latest_run_cmd(parameters, miseq_in, miseq_out)
        log.write('LATEST_RUNS\t' + ' '.join(latest_runs_cmd) + '\n')
        subprocess.call(latest_runs_cmd)
        os.remove(NON_EXCLUDED_RUNS.name)

# Computing the dump files
def compute_dump_files(parameters, log):
    MANIFESTS = list(utils.get_manifests(parameters).keys())
    # Input files indicating directory of latest run
    MISEQ_LATEST_RUNS = utils.get_MiSeq_latest_runs_file(parameters)
    # Output files
    DUMP_FILES = open(utils.get_dump_tsv_file(parameters), 'w')
    DUMP_FILES.write(DUMP_HEADER)
    # List of runs for which the dump file will be generated
    RUNS_TO_DO = []
    for manifest in MANIFESTS:
        for _, run_data in pd.read_csv(MISEQ_LATEST_RUNS[manifest], sep=',').iterrows():
            run, latest_run = run_data[RUN_ID], utils.get_latest_run_dir(run_data)
            if pd.isnull(latest_run):
                dump_cmd_out, dump_out = 'no_latest_run', str(np.nan)
            else:
                dump_out = utils.get_dump_out(parameters, run)
                if not os.path.isfile(dump_out):
                    dump_cmd = utils.get_dump_cmd(parameters, run, latest_run)
                    subprocess.call(dump_cmd)
                    dump_cmd_out = ' '.join(dump_cmd)
                else:
                    dump_cmd_out = 'file exists'
            DUMP_FILES.write('\n' + run + '\t' + dump_out)
            log.write('DUMP_FILES\t' + run + '\t' + dump_cmd_out + '\t' + dump_out + '\n')
    DUMP_FILES.close()

# Extracting matched variants into a csv file
def compute_variants(parameters, log):

    def get_run_list(runs_df):
        '''Extract all runs to consider from dataframe of runs'''
        RUN_LIST = []
        for _, run_data in runs_df.iterrows():
            RUN_LIST.append(run_data[RUN_ID])
        return RUN_LIST

    def get_excluded_samples(excluded_samples_df, runs, manifest):
        '''Get excluded samples for a manifest and a list of runs'''
        EXCLUDED_SAMPLES = {}
        for run in runs:
            # By default, no excluded sample per run
            EXCLUDED_SAMPLES[run] = []
        for _, sample_data in excluded_samples_df.iterrows():
            # Sample = integer DNA-19157-CG001Qv40Run179-4 -> 4
            run, sample = sample_data[EXC_RUN], sample_data[EXC_SAMPLE].split('-')[-1]
            if manifest in run and run in runs:
                EXCLUDED_SAMPLES[run].append(sample)
        return EXCLUDED_SAMPLES

    def get_MSI_coordinates(manifest_file):
        '''Get a dictionary indexed by chromosome of a list (start,end) of amplicons
        coordinates for amplicons labeled as MSI amplicons'''
        amplicons_df = pd.read_csv(manifest_file, sep='\t')
        MSI_df = amplicons_df.loc[amplicons_df[AMP_MUT]==AMP_MSI]
        MSI_coords = defaultdict(list)
        for _, amplicon in MSI_df.iterrows():
            MSI_coords[amplicon[AMP_CHR]].append((int(amplicon[AMP_START]), int(amplicon[AMP_END])))
        return MSI_coords

    def test_out_of_MSI(variant, MSI_coords):
        '''True if variant not in an MSI amplicon'''
        chromosome, position = variant[DUMP_CHR], variant[DUMP_POS]
        if chromosome in MSI_coords.keys():
            for (a_start, a_end) in MSI_coords[chromosome]:
                if a_start <= position <= a_end: return False
        return True

    def add_variant(variants_dict, variant,run, excluded_samples, MSI_coords):
        '''Add variant to variants_dict if it is in the target mutations (indels)
        not in an excluded sample and not in an MSI amplicon
        Returns 0 if variant added and 1 if filtered out'''
        if variant[DUMP_MUT] in DUMP_MUT_TO_SELECT:
            sample_run = variant[DUMP_SAMPLE].replace(DUMP_INDEL_SUFFIX, '')
            sample_1 = sample_run.split('-CG001')[0]
            sample_2 = sample_run.split('_S')[1]
            test_not_excluded_sample = sample_2 not in excluded_samples[run]
            test_not_in_MSI = test_out_of_MSI(variant, MSI_coords)
            if test_not_excluded_sample and test_not_in_MSI:
                sample = sample_1 + '_' + sample_2
                index = '_'.join([variant[DUMP_CHR], str(variant[DUMP_POS]), variant[DUMP_REF], variant[DUMP_ALT]])
                v_id = '_'.join([run, sample, index])
                variants_dict[v_id]['index'] = index
                variants_dict[v_id]['sample'] = sample_run
                variants_dict[v_id]['gene'] = variant[DUMP_GENE]
                variants_dict[v_id]['codon'] = variant[DUMP_CODON]
                variants_dict[v_id]['cDNA_change'] = variant[DUMP_CDNA]
                variants_dict[v_id]['MiSeq_run'] = run
                variants_dict[v_id]['sample'] = sample
                variants_dict[v_id]['vaf'] = variant[DUMP_VAF]
                variants_dict[v_id]['coverage'] = variant[DUMP_COV]
                variants_dict[v_id]['score'] = variant[DUMP_SCORE]
                return 0
            else: return 1
        else: return 0

    def add_variants_run(run_data, dump_files_dict, excluded_samples, MSI_coords, variants_dict, log):
        '''Add to variants_dict all variants for run pair run_pair'''
        miseq_run = run_data[RUN_ID]
        miseq_dump = dump_files_dict[miseq_run]
        if pd.isnull(miseq_dump):
            log.write('VARIANTS\t' + miseq_run + '\t' + 'missing dump file\n')
        else:
            nb_excluded_variants = 0
            for _, variant in pd.read_csv(miseq_dump, sep='\t').iterrows():
                nb_excluded_variants += add_variant(variants_dict, variant, miseq_run, excluded_samples,  MSI_coords)
            log.write('VARIANTS\t' + miseq_run + '\t' + str(nb_excluded_variants) + 'excluded variants\n')

    MANIFESTS_INFO = utils.get_manifests(parameters)
    DUMP_FILES_DF = pd.read_csv(utils.get_dump_tsv_file(parameters), sep='\t')
    DUMP_FILES_DICT = {x[DUMP_RUN]: x[DUMP_FILE] for _, x in DUMP_FILES_DF.iterrows()}
    MISEQ_LATEST_RUNS_FILES = utils.get_MiSeq_latest_runs_file(parameters)
    EXCLUDED_SAMPLES_FILES = utils.get_excluded_samples_file(parameters)
    VARIANTS_CSV_FILES = utils.get_variants_csv_file(parameters)
    for manifest, manifest_file in MANIFESTS_INFO.items():
        MSI_COORDS = get_MSI_coordinates(manifest_file)
        RUNS_DF = pd.read_csv(MISEQ_LATEST_RUNS_FILES[manifest], sep=',')
        RUNS_LIST = get_run_list(RUNS_DF)
        EXCLUDED_SAMPLES_DF = pd.read_csv(EXCLUDED_SAMPLES_FILES[manifest], sep=',')
        EXCLUDED_SAMPLES = get_excluded_samples(EXCLUDED_SAMPLES_DF, RUNS_LIST, manifest)
        VARIANTS_DICT = defaultdict(dict)
        for _, run_data in RUNS_DF.iterrows():
            add_variants_run(run_data, DUMP_FILES_DICT, EXCLUDED_SAMPLES, MSI_COORDS, VARIANTS_DICT, log)
        VARIANTS_DF = pd.DataFrame.from_dict(VARIANTS_DICT, orient='index')
        VARIANTS_DF.to_csv(VARIANTS_CSV_FILES[manifest], sep=',')

# Reading parameters
PARAMETERS = utils.read_parameters(sys.argv[1])
# Log file
LOG = open(sys.argv[2], 'w')
# Commands to run
CMDS = sys.argv[3:]

if 'excluded_runs' in CMDS:
    filter_excluded_runs(PARAMETERS, LOG)
if 'latest_runs' in CMDS:
    compute_latest_runs(PARAMETERS, LOG)
if 'dump_files' in CMDS:
    compute_dump_files(PARAMETERS, LOG)
if 'variants' in CMDS:
    compute_variants(PARAMETERS, LOG)

LOG.close()
