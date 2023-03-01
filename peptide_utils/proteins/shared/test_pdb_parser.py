import re, random
from pathlib import Path
import json
import pandas as pd
import numpy as np
from pandas import DataFrame
import os, sys
sys.path.append(os.path.abspath('.'))
from peptide_utils.proteins.shared.protein import pdb_to_string, from_pdb_file, get_seq_from_pdb_file_pdb_parser
from peptide_utils.proteins.shared.residue_constants import id_to_restype
from icecream import ic
ic.configureOutput(includeContext=True, argToStringFunction=lambda _: str(_))
from utils.log_util import logger
from peptide_utils.pdb_utils import PDB_ID, PEP_CHAIN, PROT_CHAIN, PEP_SEQ, PROT_SEQ
from utils.log_util import get_logger
import logging
from tqdm import tqdm


logger = get_logger(name=__name__, log_file='peptide_utils/proteins/shared/test.log', log_level=logging.INFO)
leiyipin_dir_v2 = Path('/mnt/sdb/bio_drug_corpus/ppi/ppi-leiyipin/peptide_protein_interaction_v2')
pdb_files_dir = Path('data/pdb/pdb_files')
version='2.3'
pos_only_seqs_pair_df_file = leiyipin_dir_v2 / f'peptide_protein_interaction_v{version}_unique.csv'


def teset_get_seq_from_pdb_file_pdb_parser():
    """  """
    pdb_filename = 'data/pdb/pdb_files_tmp/1en2.pdb'
    chain = 'A'
    seq = get_seq_from_pdb_file_pdb_parser(pdb_filename, chain_id=chain)
    logger.info(f'{pdb_filename} chain {chain}')
    logger.info(f'seq {seq}, len(seq) {len(seq)}')
    return seq


special_seq_pdb_id_and_chains = (
    # the ending of H chain IDQFGCSSC vs IDQFGCSSVLIVVC(right), it is wrong to split()!
    #  ATOM   1821  CD2BLEU I  16      74.274  43.196  36.635  0.20 12.89           C
    ('1h8d', 'H'),  

    # ATOM   1889  SG ACYS I  27      79.266  36.171  35.578  0.75  7.75           S  
    # ATOM   1890  N  BSER I  27      82.342  35.812  36.571  0.25 15.65           N  
    ('1h9h', 'I'),
    ('1h9h', 'E'),

)


def show_diff(seq1, seq2, row_i, pdb_id, chain):
    """  """
    if seq1 != seq2:
        logger.info(f'row_i {row_i}\n{pdb_id} {chain} len(seq1) {len(seq1)}\n{seq1}\n{seq2}')
        for char_index, (char1, char2) in enumerate(zip(seq1, seq2)):
            if char1 != char2:
                ic(f'{char_index} {char1} {char2}')
                break
        raise Exception('seq diff')


def main():
    """  
    """
    logger.info('starts')
    special_pdb_id = ''.lower()
    logger.info(f'special_pdb_id {special_pdb_id}')

    seqs_pair_df = pd.read_csv(pos_only_seqs_pair_df_file)
    ic(len(seqs_pair_df))
    ic(seqs_pair_df.columns.tolist())
    for row_i, row in tqdm(seqs_pair_df.iterrows()):
        # if row_i < 2: continue
        # if row_i > 2: break
        # if row_i == 2: ic(row)
        
        pdb_id, pep_chain, prot_chain, pep_seq, prot_seq = row.values.tolist()
        if special_pdb_id and pdb_id != special_pdb_id: continue
        
        pdb_filename = pdb_files_dir / f'{pdb_id}.pdb'
        if (pdb_id, pep_chain) not in special_seq_pdb_id_and_chains:
            pep_seq2 = get_seq_from_pdb_file_pdb_parser(pdb_filename, pep_chain)
            show_diff(pep_seq, pep_seq2, row_i, pdb_id, pep_chain)
            
        if (pdb_id, prot_chain) not in special_seq_pdb_id_and_chains:
            prot_seq2 = get_seq_from_pdb_file_pdb_parser(pdb_filename, prot_chain)
            show_diff(prot_seq, prot_seq2, row_i, pdb_id, prot_chain)

        if special_pdb_id and pdb_id == special_pdb_id: break
    logger.info('end')

# main()
teset_get_seq_from_pdb_file_pdb_parser()