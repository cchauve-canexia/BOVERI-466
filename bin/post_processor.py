import os
import re
from intervaltree import IntervalTree
from common import codon_converter
import warnings
from bunch import Bunch
BASE_DIR = os.path.dirname(os.path.realpath(__file__))
SNPEFF_PRIORITIES_PATH = os.path.join(BASE_DIR, 'snpeff_4.3_priorities.txt')
SNPEFF_PRIORITIES = {}
with open(SNPEFF_PRIORITIES_PATH, 'r') as handle:
    ind = 0
    for line in handle.readlines():
        line = line.strip()
        if len(line) > 0:
            SNPEFF_PRIORITIES[line] = ind
            ind += 1


def is_hotspot(hotspot_intervals, chrom, pos, ref, alt):
    """
    @hotspot_intervals: dict with chromosome keys and IntervalTree values
    @chrom: str
    @pos: int
    @ref: str
    @alt: str
    returns boolean whether the target overlaps with hotspot intervals
    """
    if chrom not in hotspot_intervals:
        return False
    else:
        end_pos = pos + max(len(ref), len(alt))
        return hotspot_intervals[chrom].overlaps(pos, end_pos)


def get_default_msa(annotations):
    """
    @annotations: list of Bunch
    returns the most significant annotation given a list of them
    """
    best_rank = 999
    if len(annotations) < 1:
        return None

    if len(annotations) == 1:
        return annotations[0]

    best_ann = annotations[0]
    if 'impact' in annotations[0] and\
        'annotation_type' in annotations[0]:
        for ann in annotations:
            key = ann['annotation_type'] # new version
            if ',' in SNPEFF_PRIORITIES.keys()[0]: # old format:
                key = '{},{}'.format(ann['annotation_type'], ann['impact'])

            if key is None:
                warnings.warn('Annotation key is null')
            else:
                rank = get_rank(key)
                if rank < best_rank and ann['coding']:
                    # lower rank is better
                    best_rank = rank
                    best_ann = ann
        return best_ann
    return None


def get_rank(keys):
    """ keys: a string where the keys are separated by '&'
    returns the rank as integer. Lower rank means more significant
    """
    best_rank = 999
    for key in keys.split("&"):
        if key in SNPEFF_PRIORITIES:
            rank = SNPEFF_PRIORITIES[key]
            if rank < best_rank:
                best_rank = rank
        else:
            warnings.warn('Annotation key {} not found '\
                'in priorities'.format(key))
    return best_rank


def create_snp_ref(snp_file):
    """
    returns an snp dictionary. Key is: chrom:pos:ref:altallele
    @snp_file: the path of the snp processed vcf file
    """
    snpref = {}
    for snpfileline in open(snp_file, 'r'):
        snpfiledata = re.split('\t', snpfileline.rstrip())
        if re.match("##", snpfileline):
            # info line so skip
            continue

        if snpfiledata[0] == '#CHROM':
            header = snpfiledata
            continue

        else:
            # create data dictionary using header as keys
            datadict = dict(zip(header, snpfiledata))

        chrom = datadict['#CHROM']
        pos = datadict['POS']
        ref = datadict['REF']
        alt = datadict['ALT']
        rsid = datadict['ID']
        info = datadict['INFO']
        infocaf = re.match("CAF=(\[.+\]);.+", info)

        # split the alt alleles and create individual rsid refs for each
        alt_alleles = re.split(',', alt)
        for altallele in alt_alleles:
            key = ':'.join([chrom, pos, ref, altallele])
            if infocaf:
                snpref[key] = {'rs':rsid,
                               'caf':infocaf.group(1)}
            else:
                snpref[key] = {'rs':rsid}

    return snpref


def add_snp_id(snpxref, dataset):
    """
    adds 'snp_id' to the dataset
    """
    chrom = str(dataset['CHROM'])
    pos = str(dataset['POS'])
    ref = dataset['REF']
    alt = dataset['ALT']
    key = ':'.join([chrom, pos, ref, alt])
    if key in snpxref:
        snp_id = snpxref[key]['rs']
        dataset['snp_id'] = snp_id


class SnpEff(object):
    @staticmethod
    def _curate_strelka_data(dataset):
        """
        @dataset: a dict
        This method sets some of the attributes that are needed but are not already
        in the strelka vcf file such as NA and NR which need to be deduced
        """
        for key in ['NORMAL', 'TUMOR']:
            dataset[key] = dataset[key].strip()

        format = dataset['FORMAT'].strip().split(':')
        normal = dataset['NORMAL'].strip().split(':')
        tumour = dataset['TUMOR'].strip().split(':')

        # depth = TA + TR or NA + NR, does not include TO/NO
        tumour_depth = int(tumour[format.index('DP')]) # Read depth for tier1
        normal_depth = int(normal[format.index('DP')])

        # Reads strongly supporting indel allele for tiers 1,2
        tir_ind = format.index('TIR')
        if tir_ind < 0:
            raise Exception('Could not find TIR in {}'.format(format))

        ta_t1_t2 = tumour[tir_ind]
        ta = int(ta_t1_t2.split(',')[0])
        to = int(tumour[format.index('TOR')].split(',')[0])
        tr = tumour_depth - ta - to
        dataset['TA'] = ta
        dataset['TR'] = tr
        dataset['TO'] = to

        # Reads strongly supporting indel allele for tiers 1,2:
        na_t1_t2 = normal[tir_ind]
        na = int(na_t1_t2.split(',')[0])
        no = int(normal[format.index('TOR')].split(',')[0])
        nr = normal_depth - na - no
        dataset['NA'] = na
        dataset['NR'] = nr
        dataset['NO'] = no
        
        #ta_t1_t2 = tumour[tir_ind]
        #ta = int(ta_t1_t2.split(',')[0])
        #tr = tumour_depth - ta
        #to = int(tumour[format.index('TOR')].split(',')[0])
        #dataset['TA'] = ta
        #dataset['TR'] = tr
        #dataset['TO'] = to
#
        ## Reads strongly supporting indel allele for tiers 1,2:
        #na_t1_t2 = normal[tir_ind]
        #na = int(na_t1_t2.split(',')[0])
        #nr = normal_depth - na
        #no = int(normal[format.index('TOR')].split(',')[0])
        #dataset['NA'] = na
        #dataset['NR'] = nr
        #dataset['NO'] = no

        return dataset

    @staticmethod
    def _convert_snpeff_ann(dataset, ann_fields):
        """ uses the ANN field and converts it into a list of dict
            in the dataset.annotations
        @dataset: dict
        @ann_fields: list of str
        """
        if 'ANN' in dataset:
            parts = dataset['ANN'].split(',')
            annotations = []
            for ann in parts:
                data = dict(zip(ann_fields, ann.split('|')))
                ann = {'annotation_type': data['Annotation'],\
                       'impact': data['Annotation_Impact'],\
                       'gene': data['Gene_Name'],\
                       'cdna_change': data['HGVS.c'],\
                       'protein_change': data['HGVS.p'],\
                       'exon': data['Rank'],\
                       'transcript': '',\
                       'coding': data['Transcript_BioType'] in ['Coding', 'protein_coding']}

                if data['Feature_Type'] == 'transcript':
                    ann['transcript'] = data['Feature_ID']

                # snpeff 4.1 uses 'Coding' while snpeff 4.3 uses 'protein_coding'
                pro = ann['protein_change'].replace('p.', '')
                ann['codon'] = codon_converter.shorten_codon_names(pro)
                annotations.append(ann)

            dataset['annotations'] = annotations


    def process(self, vcf, hotspot_intervals, target_interval):
        """ processes the vcf.data which is a list of dict
        """
        for dataset in vcf.data:
            self._convert_snpeff_ann(dataset, vcf.annotation_headers)
            chrom = str(dataset['CHROM'])
            if isinstance(chrom, int):
                dataset['CHROM'] = "chr%i" % chrom
            elif not chrom.startswith('chr'):
                dataset['CHROM'] = "chr%s" % chrom

            dataset['HOTSPOT'] = is_hotspot(
                hotspot_intervals, dataset['CHROM'], dataset['POS'],
                dataset['REF'], dataset['ALT'])
            if vcf.source == 'strelka':
                dataset = SnpEff._curate_strelka_data(dataset)


class TransVar(SnpEff):
    def _merge_data(self, data):
        """
        @data: list of dict
        This method merges the list of dict into a shorter list
            grouping on the chromosome and position
        returns a list of dict
        """
        results = {}
        for dataset in data:
            key = "{}:{}".format(dataset['CHROM'], dataset['POS'])
            gdna, cdna, protein = dataset['coordinates(gDNA/cDNA/protein)'].split('/')
            exon = dataset['region'].split('_exon_')[-1].replace('_and_', ',')
            tr = dataset['transcript'].split(' ')[0]
            ann = {'gene': dataset['gene'],\
                   'transcript': tr,\
                   'cdna_change': cdna,\
                   'protein_change': protein,\
                   'exon': exon\
                  }

            if key in results:
                results[key]['annotations'].append(ann)
            else:
                results[key] = dataset
                results[key]['annotations'] = [ann]

        return results.values()


    def process(self, vcf, hotspots, target_interval):
        """ processes the vcf.data which is a list of dict
        """
        vcf.data = self._merge_data(vcf.data)
        super(TransVar, self).process(vcf, hotspots, target_interval)


def get_msa(annotations, tr_ids):
    """
    @annotations: dict
    @tr_ids list of transcript ids without the version
    """
    transcript_annotations = [ann for ann in annotations if ann['transcript'].split('.')[0] in tr_ids ]
    if transcript_annotations:
        return get_default_msa(transcript_annotations)
    return get_default_msa(annotations)


def get_hotspot_intervals(hotspots):
    """
    returns a dict with chromosomes as keys and IntervalTree() values
    :param hotspots: a dict returned from the function: manifest_controller.get_hotspots() that contains amplicon ids as
        keys, and Bunch values. The Bunch is constructed from data found in the hotspot manifest file.
    """
    intervals = {} # keys are chromosomes, values are Intervaltree()
    for obj in hotspots.values():
        if not obj.chrom in intervals:
            intervals[obj.chrom] = IntervalTree()
        # add 1 to end to make it inclusive
        intervals[obj.chrom][obj.start:obj.end + 1] = obj
    return intervals


def __add_more_fields(msa, dataset, source):
    msa.chromosome = dataset.get('CHROM', '')
    msa.pos = int(dataset.get('POS', 0))
    msa.ref = dataset.get('REF', '')
    msa.alt = dataset.get('ALT', '')
    msa.coverage = dataset.get('TA', 0) + dataset.get('TR', 0) + dataset.get('TO', 0)
    if msa.coverage != 0:
        msa.vaf = float(dataset.get('TA', 0.0)) / msa.coverage

    msa.mutation_type = 'snv'
    if 'strelka' in source.lower():
        msa.mutation_type = 'indel'

    msa.normal_ref_allele_count = dataset.get('NR', 0)
    msa.normal_alt_allele_count = dataset.get('NA', 0)
    msa.normal_oth_allele_count = dataset.get('NO', 0)
    msa.tumour_ref_allele_count = dataset.get('TR', 0)
    msa.tumour_alt_allele_count = dataset.get('TA', 0)
    msa.tumour_oth_allele_count = dataset.get('TO', 0)
    msa.score = dataset.get('PR', 0.0) # Probability of somatic mutation
    if 'QSI' in dataset:
        msa.score = dataset['QSI']
    msa.filter = dataset.get('FILTER', '')


def process_annotations(vcf, hotspots, annotator_type='snpeff'):
    """ processes the vcf file which is
    adding additional fields to dicts inside vcf.data
    :param vcf: Vcf object
    :param hotspots: a dict returned from the function: manifest_controller.get_hotspots() that contains amplicon ids as
        keys, and Bunch values. The Bunch is constructed from data found in the hotspot manifest file.
    :param annotator_type: A string that determines which annotator to use. One of 'snpeff' or 'transvar'.
    """
    annotator = None
    if annotator_type == 'transvar':
        annotator = TransVar()
    elif annotator_type == 'snpeff':
        annotator = SnpEff()
    else:
        raise Exception('Unsupported annotator: {}'.format(annotator_type))

    transcript_ids = [val.transcript_id.split(".")[0] for val in hotspots.values()]
    hotspot_intervals = get_hotspot_intervals(hotspots)
    annotator.process(vcf, hotspot_intervals, hotspot_intervals)
    for val in vcf.data:
        msa = Bunch(get_msa(val['annotations'], transcript_ids))
        __add_more_fields(msa, val, vcf.source)
        val['msa'] = msa


