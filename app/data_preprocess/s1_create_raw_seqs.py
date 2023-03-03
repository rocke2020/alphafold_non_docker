import re, random
from pathlib import Path
import json
import pandas as pd
import numpy as np
from pandas import DataFrame
import os, sys
from icecream import ic
ic.configureOutput(includeContext=True, argToStringFunction=lambda _: str(_))
sys.path.append(os.path.abspath('.'))
from peptide_utils.pdb_utils import get_pdb_seq_from_pdb_file, is_hunman_type_by_checking_pdb
from utils.log_util import logger
from peptide_utils.pep_utils import is_natural_only_supper
from data_process.datasets.uniprot.common import get_uniprot_seq_from_pdb_id_and_chain


test_data_dir = Path('/mnt/sda/bio_drug_corpus/ppi/Improving Peptide-Protein Docking with AlphaFold-Multimer using Forced Sampling')
gold_positive_test_file = test_data_dir / 'test-complex-num112-all positive.CSV'
gold_negative_test_file = test_data_dir / 'test-complex-num672-1pos5neg.CSV'

pdb_seqs_dir = test_data_dir / 'pdb_seqs'
pdb_seqs_dir.mkdir(exist_ok=1)
gold_pos_pdb_seqs_file = pdb_seqs_dir / 'positive112.csv'
gold_neg_pdb_seqs_file = pdb_seqs_dir / 'negative560.csv'
gold_pos_pdb_seqs_natural_aa_file = pdb_seqs_dir / 'positive112_natural_aa.csv'
gold_neg_pdb_seqs_natural_aa_file = pdb_seqs_dir / 'negative560_natural_aa.csv'

fasta_seqs_dir = test_data_dir / 'fasta_seqs'
fasta_seqs_dir.mkdir(exist_ok=1)
gold_pos_pdb_fasta_seqs_file = fasta_seqs_dir / 'positive112.csv'
gold_neg_pdb_fasta_seqs_file = fasta_seqs_dir / 'negative560.csv'
gold_pos_pdb_fasta_seqs_natural_aa_file = fasta_seqs_dir / 'positive112_natural_aa.csv'
gold_human_pos_pdb_fasta_seqs_natural_aa_file = fasta_seqs_dir / 'positive112_natural_aa_human.csv'
gold_human_mid_len_pos_pdb_fasta_seqs_natural_aa_file = fasta_seqs_dir / 'positive112_natural_aa_human_200_300mer.csv'
gold_neg_pdb_fasta_seqs_natural_aa_file = fasta_seqs_dir / 'negative560_natural_aa.csv'

uniprot_seqs_dir = test_data_dir / 'uniprot_seqs'
uniprot_seqs_dir.mkdir(exist_ok=1)

pdb_fasta_seqres_dict_file = Path('data/pdb/pdb_seqres.json')
with open(pdb_fasta_seqres_dict_file, 'r', encoding='utf-8') as f:
    pdb_fasta_seqres_dict = json.load(f)

version='2.3'
leiyipin_dir = Path('/mnt/sda/bio_drug_corpus/ppi/ppi-leiyipin/peptide_protein_interaction_v2')
pos_only_seqs_pair_df_file = leiyipin_dir / f'peptide_protein_interaction_v{version}_unique.csv'
camp_pos_pdb_fasta_seqs_file = fasta_seqs_dir / f'camp_pos_pdb_fasta_seqs.csv'
camp_pos_pdb_fasta_seqs_uniprot_file = fasta_seqs_dir / f'camp_pos_pdb_fasta_seqs_uniprot.csv'


def create_pdb_seqs(df:DataFrame, out_file:Path, select_positive=True):
    """  """
    prot_seqs, pep_seqs = create_seqs(df, select_positive, get_pdb_seq_from_pdb_file)
    df['pep_pdb_seq'] = pep_seqs
    df['prot_pdb_seq'] = prot_seqs
    df.to_csv(out_file, index=False, sep=',')

    natural_aa_df = df[
        (df['pep_pdb_seq'].map(is_natural_only_supper)) & (df['prot_pdb_seq'].map(is_natural_only_supper))]
    ic(len(natural_aa_df))
    _out_file = out_file.with_stem(f'{out_file.stem}_natural_aa')
    natural_aa_df.to_csv(_out_file, index=False, sep=',')


def create_pdb_fasta_seqs(df:DataFrame, out_file:Path, select_positive=True):
    """  """
    prot_seqs, pep_seqs = create_seqs(df, select_positive, get_pdb_fasta_seq_from_dict)
    df['pep_pdb_fasta_seq'] = pep_seqs
    df['prot_pdb_fasta_seq'] = prot_seqs
    ic(len(df))
    df.to_csv(out_file, index=False, sep=',')

    natural_aa_df = df[
        (df['pep_pdb_fasta_seq'].map(is_natural_only_supper)) & (df['prot_pdb_fasta_seq'].map(is_natural_only_supper))]
    ic(len(natural_aa_df))
    _out_file = out_file.with_stem(f'{out_file.stem}_natural_aa')
    natural_aa_df.to_csv(_out_file, index=False, sep=',')


def create_seqs(df, select_positive, get_seq_func):
    prot_seqs = []
    pep_seqs = []
    for i, row in df.iterrows():
        prot_chain = row['prot_chain']
        pep_chain = row['pep_chain']
        if select_positive:
            pdb_id_prot = row['pdb_id']
            pdb_id_pep = pdb_id_prot
        else:
            pdb_id_prot = row['pdb_id_prot']
            pdb_id_pep = row['pdb_id_pep']
        
        prot_seq = get_seq_func(pdb_id_prot, prot_chain)
        prot_seqs.append(prot_seq)
        pep_seq = get_seq_func(pdb_id_pep, pep_chain)
        pep_seqs.append(pep_seq)
    return prot_seqs, pep_seqs


def get_pdb_fasta_seq_from_dict(pdb_id, chain):
    """  """
    pdb_id_chain = pdb_id.lower() + '_' + chain
    seq  = pdb_fasta_seqres_dict.get(pdb_id_chain, '')
    if seq == '':
        logger.info(f'pdb_id_chain {pdb_id_chain} not fasta seq')
    return seq


def create_seqs_gold_test112():
    """ 112 positive and 560 negative """
    # 'pdb_id_prot', 'prot_chain', 'pdb_id_pep', 'pep_chain'
    gold_neg_test_df = pd.read_csv(gold_negative_test_file)
    ic(gold_neg_test_df.columns)
    ic(gold_neg_test_df.shape)
    
    # 'pdb_id', 'prot_chain', 'pep_chain'
    gold_pos_test_df = pd.read_csv(gold_positive_test_file)
    ic(gold_pos_test_df.columns)
    ic(gold_pos_test_df.shape)

    create_pdb_fasta_seqs(gold_neg_test_df, gold_neg_pdb_fasta_seqs_file, select_positive=False)
    create_pdb_fasta_seqs(gold_pos_test_df, gold_pos_pdb_fasta_seqs_file)
    logger.info('Finishes to create pdb_fasta_seqs')

    create_pdb_seqs(gold_neg_test_df, gold_neg_pdb_seqs_file, select_positive=False)
    create_pdb_seqs(gold_pos_test_df, gold_pos_pdb_seqs_file)


def create_fasta_seqs_for_camp_positive():
    """  """
    unique_positive_data = pd.read_csv(pos_only_seqs_pair_df_file)
    ic(unique_positive_data.columns)
    unique_positive_data.rename(
        columns={'Peptide Chain': 'pep_chain', 'Receptor Chain': 'prot_chain', 'PDB': 'pdb_id'}, inplace=True)
    create_pdb_fasta_seqs(unique_positive_data, camp_pos_pdb_fasta_seqs_file)
    logger.info('Finishes to create pdb_fasta_seqs')


def create_uniprot_seqs(input_file: Path):
    """ input_file is a saved pandas dataframe file """
    ic(input_file.name)
    df = pd.read_csv(input_file)
    ic(df.columns)
    ic(len(df))
    df['uniprot_seq'] = df.apply(create_unpprot_seq_for_prot_chain, axis=1)
    df = df[df['uniprot_seq'].map(lambda x: len(x) > 0)]
    ic(len(df))
    df = df[df['uniprot_seq'].map(is_natural_only_supper)]
    ic(len(df))
    ic(df.head())
    df['uniprot_seq_len'] = df['uniprot_seq'].map(len)
    ic(df['uniprot_seq_len'].tolist())
    max_uniprot_seq_len = max(df['uniprot_seq_len'])
    min_uniprot_seq_len = min(df['uniprot_seq_len'])
    ic(max_uniprot_seq_len)
    ic(min_uniprot_seq_len)
    # py38 has no with_stem, but py39 has
    out_file = input_file.with_name(f'{input_file.stem}_uniprot.csv')
    df.to_csv(out_file, index=False, sep=',')


def create_unpprot_seq_for_prot_chain(row):
    prot_chain = row['prot_chain']
    if 'pdb_id' in row:
        pdb_id_prot = row['pdb_id']
    else:
        pdb_id_prot = row['pdb_id_prot']
    seq = get_uniprot_seq_from_pdb_id_and_chain(pdb_id_prot, prot_chain)
    return seq
    

def check_fasta_seqs():
    """  """
    pos_pdb_fasta_seq_df = pd.read_csv(gold_pos_pdb_fasta_seqs_file)
    pos_pdb_fasta_seq_natural_df = pd.read_csv(gold_pos_pdb_fasta_seqs_natural_aa_file)
    neg_pdb_fasta_seq_df = pd.read_csv(gold_neg_pdb_fasta_seqs_file)
    neg_pdb_fasta_seq_natural_df = pd.read_csv(gold_neg_pdb_fasta_seqs_natural_aa_file)

    ic(neg_pdb_fasta_seq_natural_df.columns)


def select_gold_human_pdb_fasta_seqs():
    """  """
    gold_pos_pdb_fasta_seqs_df = pd.read_csv(gold_pos_pdb_fasta_seqs_natural_aa_file)
    ic(gold_pos_pdb_fasta_seqs_df.columns)
    ic(len(gold_pos_pdb_fasta_seqs_df))
    human_df = gold_pos_pdb_fasta_seqs_df[
        gold_pos_pdb_fasta_seqs_df['pdb_id'].map(is_hunman_type_by_checking_pdb)].copy()
    ic(len(human_df))
    human_df.to_csv(gold_human_pos_pdb_fasta_seqs_natural_aa_file, index=False, sep=',')

    human_df_mid_len_df = human_df[
        (human_df['prot_pdb_fasta_seq'].map(lambda x: len(x) >=200)) & 
        (human_df['prot_pdb_fasta_seq'].map(lambda x: len(x) <=300))]
    ic(len(human_df_mid_len_df))
    human_df_mid_len_df.to_csv(gold_human_mid_len_pos_pdb_fasta_seqs_natural_aa_file, index=False, sep=',')

    
if __name__ == "__main__":
    # create_seqs_gold_test112()
    # check_fasta_seqs()
    # create_fasta_seqs_for_camp_positive()
    # select_gold_human_pdb_fasta_seqs()
    create_uniprot_seqs(gold_human_mid_len_pos_pdb_fasta_seqs_natural_aa_file)
    ic('end')