from collections import defaultdict
from pandas import Series, DataFrame
import random, re
import numpy as np
from icecream import ic
ic.configureOutput(includeContext=True)
from pathlib import Path
from utils.log_util import get_logger
import logging
from tqdm import tqdm


logger = get_logger(name=__name__, log_file=None, log_level=logging.INFO)
SEQUENCE = 'Sequence'
# 20 natrual
aminoacids = ["A", "C", "D", "E", "F", "G", "H", "I", "L", "M", "N", "P", "K", "Q", "R", "S", "T", "V", "W", "Y"]
aminoacids_plus_BX = aminoacids + ["B", "X"]
### Non-containing methionine sequences were preferred. Methionine，简写M，Met, 甲硫氨酸, for diversity, not used in codes.
# often_AMP_aminoacids = aminoacids.copy()
# often_AMP_aminoacids.remove('M')
C_TERMINUS = "cTerminus"
N_TERMINUS = "nTerminus"
LEAST_SEQ_LENGTH = 3
MAX_SEQ_LEN = 50
NOT_HEMOLYTIC_KEY = 'isNotHemolytic'
ACTIVITY_KEY = 'activity'

not_digits_at_head_end_pat = re.compile(r'^\D*|\D*$')
valid_units = ('µg/ml', 'µM')
NAs = ["NA", "na", "Na", "nA", "N/A", "n/a"]
NA = 'NA'
full_replace_pat = re.compile(r'\s+|[>=<]')
non_az_at_head_pat = re.compile(r'^[^a-zA-Z]+')
hemolysis_names = ['hemolisis', 'hemolysis']
cytotoxicity_names = ['cytotoxicity', 'cell death']


def is_seq_len_valid(seq:str):
    if LEAST_SEQ_LENGTH <= len(seq) <= MAX_SEQ_LEN:
        return True
    return False


def calc_natural_aa_ratio(seq):
    len_seq = len(seq)
    cnt = 0
    for i in seq:
        if i in aminoacids :
            cnt = cnt+1
    score = cnt / len_seq
    return score


def is_80percent_natural(seq):
    """  """
    count = 0
    for char in seq:
        if char in aminoacids:
            count += 1
    if count / len(seq) >= 0.8:
        return True
    return False


def is_natural_and_variable_length(seq, least_len=5, max_len=15):
    if seq and len(seq) >= least_len and len(seq) <= max_len and is_natural_only_supper(seq):
        return True
    return False


def is_natural_and_length_5_15(seq):
    if seq and len(seq) >= 5 and len(seq) <= 15 and is_natural_only_supper(seq):
        return True
    return False


def is_natural_and_max_length_30(seq):
    if seq and len(seq) >= LEAST_SEQ_LENGTH and len(seq) <= 30 and is_natural_only_supper(seq):
        return True
    return False

def is_natural_and_valid_length(seq):
    if seq and len(seq) >= LEAST_SEQ_LENGTH and len(seq) <= MAX_SEQ_LEN and is_natural_only_supper(seq):
        return True
    return False


def is_natural_and_valid_length_shortest4(seq):
    if seq and len(seq) >= 4 and len(seq) <= MAX_SEQ_LEN and is_natural_only_supper(seq):
        return True
    return False


def is_natural_only_supper(seq):
    """ If the char is not in upper case, treat as not natural """
    if isinstance(seq, str) and seq:
        for aa in seq:
            if aa not in aminoacids:
                return False
        return True
    return False


def get_not_natural_aas_only_supper(seq):
    """ If the char is not in upper case, treat as not natural """
    not_natural_aas = []
    if isinstance(seq, str) and seq:
        for aa in seq:
            if aa not in aminoacids:
                not_natural_aas.append(aa)
    return not_natural_aas


def randomChoice(l):
    return l[random.randint(0, len(l) - 1)]


def novelty(seqs, list_):
    novel_seq = []
    for s in seqs:
        if s not in list_:
            novel_seq.append(s)
    return novel_seq, (len(novel_seq)/len(seqs))*100


def uniqueness(seqs):
    """ Use dict to store the repetition num which is used """
    unique_seqs_dict = defaultdict(int)
    for s in seqs:
        if len(s) > 0:
            unique_seqs_dict[s] += 1
    return unique_seqs_dict, (len(unique_seqs_dict)/len(seqs))*100


def float_ignore_plus_minus(mynumber):
    try:
        return sum(map(float,mynumber.split("±")))
    except:
        return float("inf")


def check_active(unit, concentration):
    if ((unit == "µM" and concentration < 10) or (unit == "nM" and concentration < 10000)
        or (unit == "µg/ml" and concentration < 32)):
        return True


def check_inactive(unit, concentration):
    if ((unit == "µM" and concentration > 10) or (unit == "nM" and concentration > 10000)
        or (unit == "µg/ml" and concentration > 32)):
        return True


def get_terminus_names(dbaasp_peptide):
    n_terminus = dbaasp_peptide.get(N_TERMINUS)
    if n_terminus:
        n_name = n_terminus['name']
    else:
        n_name = 'nan'
    c_terminus = dbaasp_peptide.get(C_TERMINUS)
    if c_terminus:
        c_name = c_terminus['name']
    else:
        c_name = 'nan'
    return n_name, c_name


def is_valid_terminus(n_name, c_name):
    if (n_name == 'nan' or n_name == 'ACT') and (c_name == 'nan' or c_name == 'AMD'):
        return True


def sum_plus_minus(s):
    return sum(map(float, s.split("±")))


def um_to_ug_ml(concentration:float, mw):
    """ µM unit: µmole / liter;
        mw, molecule weight unit: g/mol

        Returns:
            µg/ml
    """
    return concentration * mw / 1000


def ug_ml_to_um(concentration:float, mw):
    """
    Args:
        µg/ml

        µM unit: µmole / liter
        mw, molecule weight unit: g/mol

        Returns:
            µmole / liter
    """
    return concentration * 1000 / mw


def mg_ml_to_um(concentration:float, mw):
    return concentration * 1000_000 / mw


def is_toxic_dbaasp(toxic_items, threshold=200):
    """
    toxic_items: [[concentration, unit, mw, target_cell, hemolysis_ratio, cell_death_ratio],]
    1) Any hemolytic ratio >= 20, treat hemolytic, that's toxic
    2) ALL hemolytic/cytotoxic activities against reported target species less+equal than 200 µg/ml as toxic
    """
    if len(toxic_items)==0:
        return 0

    ## All hemolytic
    # hemolytic = True
    # hemolytic_count = 0
    # for data in toxic_items:
    #     concentration, unit, mw, hemolysis_ratio = data
    #     if hemolysis_ratio != NA:
    #         hemolytic_count += 1
    #         if hemolysis_ratio < 20:
    #             hemolytic = False
    #             break
    # if hemolytic_count and hemolytic:
    #     return 1

    # Any hemolytic
    for data in toxic_items:
        concentration, unit, mw, hemolysis_ratio = data
        if hemolysis_ratio != NA and hemolysis_ratio >= 20:
            return 1

    for data in toxic_items:
        concentration = data[0]
        unit = data[1]
        mw = data[2]
        if unit == "µM":
            concentration = um_to_ug_ml(concentration, mw)
        if concentration > threshold:
            return 0
    return 1


def is_nontoxic_dbaasp(toxic_items, threshold=250):
    """ hemolytic/cytotoxic activities against ALL reported target species larger+equal than 250 µg/ml as non toxic """
    if len(toxic_items)==0:
        return 0
    for data in toxic_items:
        concentration = data[0]
        unit = data[1]
        mw = data[2]
        if unit == "µM":
            concentration = um_to_ug_ml(concentration, mw)
        if concentration < threshold:
            return 0
    return 1


def filter_short_amino_acids(df, shortest_length=5):
    df = df[df.Sequence.map(lambda x: len(x) >= shortest_length)]
    return df
