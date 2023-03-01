from pathlib import Path
from collections import defaultdict
from utils.log_util import logger
import pandas as pd
from .pep_utils import is_natural_only_supper, SEQUENCE
from rapidfuzz.string_metric import levenshtein as lev_dist
from pandas import DataFrame
from Bio import Align


uniprot_sprot_file = '/mnt/sda/bio_drug_corpus/uniport/uniprot_sprot.fasta'
uniprot_sprot_natural_file = Path('/mnt/sda/bio_drug_corpus/uniport/uniprot_sprot_natural.csv')


def get_len_seq_dict(sequences):
    """  """
    length_seq = defaultdict(list)
    for seq in sequences:
        length_seq[len(seq)].append(seq)
    return length_seq


def get_uniprot(read_cache=1):
    """ need about 6 seconds to load saved csv file
    uniprot_sprot.fasta
    >sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-001R PE=4 SV=1
    MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPS
    EKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLD
    AKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHL
    EKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDD
    SFRKIYTDLGWKFTPL
    >sp|Q6GZX3|002L_FRG3G Uncharacterized protein 002L OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-002L PE=4 SV=1
    MSIIGATRLQNDKSDTYSAGPCYAGGCSAFTPRGTCGKDWDLGEQTCASGFCTSQPLCAR
    IKKTQVCGLRYSSKGKDPLVSAEWDSRGAPYVRCTYDADLIDTQAQVDQFVSMFGESPSL
    AERYCMRGVKNTAGELVSRVSSDADPAGGWCRKWYSAHRGPDQDAALGSFCIKNPGAADC
    KCINRASDPVYQKVKTLHAYPDQCWYVPCAADVGELKMGTQRDTPTNCPTQVCQIVFNML
    DDGSVTMDDVKNTINCDFSKYVPPPPPPKPTPPTPPTPPTPPTPPTPPTPPTPRPVHNRK
    VMFFVAGAVLVAILISTVRW
    """
    if uniprot_sprot_natural_file.is_fifo() and read_cache:
        swissprot = pd.read_csv(uniprot_sprot_natural_file)
        logger.info('Loads uniprot sequences')
    else:
        seq = ''
        newvalues_dictionary = {}
        first = True
        with open(uniprot_sprot_file) as inFile:
            for line in inFile:
                line = line.strip()
                if line[0] == ">":
                    if first ==  True:
                        cid = line.replace(">", "")
                        first = False
                        continue
                    newvalues_dictionary[cid] = seq
                    cid = line.replace(">", "")
                    seq = ""
                else:
                    seq+=line
        if seq:
            newvalues_dictionary[cid] = seq
        swissprot = pd.DataFrame(newvalues_dictionary.items(), columns=['ID', 'Sequence'])
        # print(f'orig len(swissprot) {len(swissprot)}')  # 568363
        swissprot = swissprot[swissprot["Sequence"].map(is_natural_only_supper)].reset_index(drop = True)
        swissprot['length'] = swissprot['Sequence'].map(len)
        # print(f'natural_only len(swissprot) {len(swissprot)}')  # 565701
        swissprot.to_csv(uniprot_sprot_natural_file, index=False)
        logger.info(f'Saves swissprot natural at {uniprot_sprot_natural_file}')
    logger.info(f'len(swissprot) {len(swissprot) }')
    return swissprot


def read_short_uniprot_data(short_seq_fasta_file):
    """ Only short sequence 4-50 and so only one line seq """
    seq_lst = []
    with open(short_seq_fasta_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'): continue
            seq_lst.append(line.strip())
    logger.info(f'len(seq_lst) {len(seq_lst)}')
    return seq_lst


def write_seqs_to_fasta(seqs, fasta_file):
    """  """
    with open(fasta_file, 'w', encoding='utf-8') as f:
        for i, seq in enumerate(seqs):
            f.write(f'>{i}\n{seq}\n')


def create_df_from_seqs(seqs, save_file=None):
    """ df columns: ID, Sequence, length """
    # seqs may be set()
    seqs = list(seqs)
    df = pd.DataFrame({
        'ID': range(len(seqs)),
        'Sequence': seqs,
        'length': [len(seq) for seq in seqs],
    })
    logger.info(f'df head {df.describe()}')
    if save_file:
        df.to_csv(save_file, index=False)
    return df


def cal_seq_similarities_between_df(query_df:DataFrame, target_df:DataFrame):
    """ only adds new columns into query df
    Returns:
        By default, only add similarities into query_df
        query_df, target_df
    """
    target_seqs = target_df[SEQUENCE]
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    def cal_least_lev_dist(querty_seq):
        """  """
        distances = [[lev_dist(querty_seq, seq), seq] for seq in target_seqs]
        distances.sort()
        return distances[0]

    def cal_highest_align_score(querty_seq):
        """  """
        distances = [[aligner.score(querty_seq, seq), seq] for seq in target_seqs]
        distances.sort(reverse=True)
        return distances[0]

    query_df['lev_dist_and_seq'] = query_df[SEQUENCE].map(cal_least_lev_dist)
    query_df['least_lev_dist_seq'] = query_df['lev_dist_and_seq'].map(lambda x: x[1])
    query_df['length_of_least_dist'] = query_df['least_lev_dist_seq'].map(len)
    query_df['value_of_least_dist'] = query_df['lev_dist_and_seq'].map(lambda x: x[0])

    query_df['highest_align'] = query_df[SEQUENCE].map(cal_highest_align_score)
    query_df['highest_align_seq'] = query_df['highest_align'].map(lambda x: x[1])
    query_df['length_of_highest_align'] = query_df['highest_align_seq'].map(len)
    query_df['value_of_highest_align'] = query_df['highest_align'].map(lambda x: x[0])

    del query_df['lev_dist_and_seq']
    del query_df['highest_align']
    return query_df
