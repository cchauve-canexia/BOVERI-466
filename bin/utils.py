'''
Auxiliary functions to manipulate parameters
'''
# Third-party imports
import pandas as pd

# Local imports
from data_format import LATEST_RUN, PIPELINE_OUTPUT

def read_parameters(parameters_file):
    '''Returns a dictionary of parameters values'''
    PARAMETERS_FILE = open(parameters_file, 'r').readlines()
    PARAMETERS = {}
    for l in PARAMETERS_FILE:
        if l[0] != '#' and len(l) > 2:
            l1 = l.rstrip().split('\t')
            PARAMETERS[l1[0]] = l1[1]
    return PARAMETERS

'''Functions to create file paths'''

def pad_file(parameters, file_key, suffix, dir_key=None):
    '''Pad a file name with a directory and a suffix'''
    if dir_key is None: dir_str = ''
    else: dir_str = parameters[dir_key] + '/'
    return dir_str + parameters[file_key] + suffix

def files_dict(parameters, file_key, suffix, dir_key=None):
    '''Returns a dictionary indexed by manifests of paths to files prefixed by a
    directory and suffixed by _manifest_suffix'''
    RESULT = {}
    for manifest in get_manifests(parameters):
        suff_str = '_' + manifest + suffix
        RESULT[manifest] = pad_file(parameters, file_key, suff_str, dir_key=dir_key)
    return RESULT

'''Manifests'''
def get_manifests(parameters):
    '''Returns dictionary manifest_id -> manifest_file'''
    DATA_DIR = parameters['DATA_DIR']
    MANIFESTS = parameters['MANIFESTS'].split(' ')
    RESULT = {}
    i = 0
    while i < len(MANIFESTS):
        RESULT[MANIFESTS[i]] = DATA_DIR + '/' + MANIFESTS[(i + 1)]
        i += 2
    return RESULT

'''Functions to access data files'''
def get_MiSeq_runs_file(parameters):
    '''Returns a dictionary indexed by manifests of paths to files containing
    the run info for MiSeq samples'''
    return files_dict(parameters, 'MISEQ_RUNS_DATA_CSV_PREFIX', '.csv', dir_key='DATA_DIR')

def get_MiSeq_excluded_runs_file(parameters):
    '''Returns a dictionary indexed by manifests of paths to files containing
    the run info for MiSeq samples'''
    return files_dict(parameters, 'EXCLUDED_RUNS_PREFIX', '.csv', dir_key='DATA_DIR')

'''Functions to create and access latest run information'''
def get_MiSeq_latest_runs_file(parameters):
    '''Returns a dictionary indexed by manifests of paths to files containing
    the latest run info for MiSeq samples'''
    return files_dict(parameters, 'MISEQ_RUNS_DATA_CSV_PREFIX', '_latest_run.csv', dir_key='OUT_DIR')

def get_latest_run_dir(run_data):
    '''Returns the path to latest run directory for a run described in the
    run_data dictionary; returns Nan if latest run information not available'''
    latest_run = run_data[LATEST_RUN]
    if pd.isnull(latest_run):
        latest_run = run_data[PIPELINE_OUTPUT]
    return latest_run

def get_latest_run_cmd(parameters, file_in, file_out):
    '''Command for running the get_pipeline_latest_output.py script'''
    RESULT = parameters['LATEST_RUNS_CMD'].replace('RUNS_IN', file_in)
    RESULT = RESULT.replace('RUNS_OUT', file_out).split(' ')
    return RESULT

'''Functions to create and access dump files'''

def get_dump_tsv_file(parameters):
    '''Returns the path to the tsv file containing the list of dump files'''
    return pad_file(parameters, 'DUMP_FILES_TSV', '', dir_key='OUT_DIR')

def get_dump_in(parameters, latest_run):
    '''Returns the -i parameters for the dump_variants script for a specific run'''
    return parameters['DUMP_IN'].replace('LATEST_RUN', latest_run)

def get_dump_out(parameters, run):
    '''Returns the -o parameters for the dump_variants script for a specific run'''
    return parameters['DUMP_OUT'].replace('RUN_NAME', run)

def get_dump_cmd(parameters, run, latest_run):
    '''Returns the command line for the dump_variants script for a specific run
    :param: run (str): run name
    :param: latest_run (str): directory to the latest run directory
    '''
    DUMP_CMD = parameters['DUMP_CMD'].split(' ')
    DUMP_IN = get_dump_in(parameters, latest_run)
    DUMP_OUT = get_dump_out(parameters, run)
    DUMP_PARAMETERS = parameters['DUMP_PARAMETERS'].split(' ')
    return DUMP_CMD + ['-d', DUMP_IN, '-o', DUMP_OUT] + DUMP_PARAMETERS

'''Functions to create the variants files'''

def get_excluded_runs_file(parameters):
    '''Returns a dictionary indexed by manifests of paths to files containing the
    list of excluded runs'''
    return files_dict(parameters, 'EXCLUDED_RUNS_PREFIX', '.csv', dir_key='DATA_DIR')

def get_excluded_samples_file(parameters):
    '''Returns a dictionary indexed by manifests of paths to files containing the
    list of excluded samples'''
    return files_dict(parameters, 'EXCLUDED_SAMPLES_PREFIX', '.csv', dir_key='DATA_DIR')

def get_variants_csv_file(parameters):
    '''Returns a dictionary indexed by manifests of paths to files containing the
    variants files for all  non-excluded runs for this manifest'''
    return files_dict(parameters, 'VARIANTS_CSV_PREFIX', '.csv', dir_key='OUT_DIR')
