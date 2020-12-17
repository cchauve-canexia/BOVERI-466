''' Fields of data files '''

# Runs info file
LATEST_RUN = 'latest_run'
PIPELINE_OUTPUT = 'Pipeline Output in S3 bucket rd-output'
RUN_ID = 'Run#'

# Dump files list and files
DUMP_RUN = 'run'
DUMP_FILE = 'dump_file'
DUMP_HEADER = '\t'.join([DUMP_RUN, DUMP_FILE])
DUMP_CHR = 'chromosome'
DUMP_POS = 'pos'
DUMP_MUT = 'mutation_type'
DUMP_INDEL = 'indel'
DUMP_MUT_TO_SELECT = [DUMP_INDEL]
DUMP_SAMPLE = 'sample'
DUMP_INDEL_SUFFIX = '_somatic_indels.vcf.snpeff.vcf'
DUMP_REF = 'ref'
DUMP_ALT = 'alt'
DUMP_GENE = 'gene'
DUMP_CODON = 'codon'
DUMP_CDNA = 'cdna_change'
DUMP_VAF = 'vaf'
DUMP_SCORE = 'score'
DUMP_COV = 'coverage'

# Excluded excluded
EXC_RUN = 'run name'
EXC_SAMPLE = 'sample name'

# Amplicons manifests
AMP_MUT = 'Mutations'
AMP_MSI = 'MSI'
AMP_CHR = 'Chr'
AMP_START = 'Start'
AMP_END = 'End'
