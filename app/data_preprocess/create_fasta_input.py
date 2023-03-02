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
from peptide_utils.pep_utils import get_not_natural_aas_only_supper
from peptide_utils.pdb_utils import get_pdb_seq_from_pdb_file
from utils.log_util import logger


version='2.3'
leiyipin_dir = Path('/mnt/sda/bio_drug_corpus/ppi/ppi-leiyipin/peptide_protein_interaction_v2')
pos_only_seqs_pair_df_file = leiyipin_dir / f'peptide_protein_interaction_v{version}_unique.csv'

out_dir = Path('/mnt/sdc/af_input')
out_dir.mkdir(exist_ok=1)
test_data_dir = Path('/mnt/sda/bio_drug_corpus/ppi/Improving Peptide-Protein Docking with AlphaFold-Multimer using Forced Sampling')
gold_positive_test_file = test_data_dir / 'test-complex-num112-all positive.CSV'
gold_negative_test_file = test_data_dir / 'test-complex-num672-1pos5neg.CSV'

test_data_dir = Path('/mnt/sda/bio_drug_corpus/ppi/Improving Peptide-Protein Docking with AlphaFold-Multimer using Forced Sampling')
pdb_seqs_dir = test_data_dir / 'pdb_seqs'
gold_pos_pdb_seqs_file = pdb_seqs_dir / 'positive112.csv'
gold_neg_pdb_seqs_file = pdb_seqs_dir / 'negative560.csv'
gold_pos_pdb_seqs_natural_aa_file = pdb_seqs_dir / 'positive112_natural_aa.csv'
gold_neg_pdb_seqs_natural_aa_file = pdb_seqs_dir / 'negative560_natural_aa.csv'

fasta_seqs_dir = test_data_dir / 'fasta_seqs'
gold_pos_pdb_fasta_seqs_file = fasta_seqs_dir / 'positive112.csv'
gold_neg_pdb_fasta_seqs_file = fasta_seqs_dir / 'negative560.csv'
gold_pos_pdb_fasta_seqs_natural_aa_file = fasta_seqs_dir / 'positive112_natural_aa.csv'
gold_neg_pdb_fasta_seqs_natural_aa_file = fasta_seqs_dir / 'negative560_natural_aa.csv'
camp_pos_pdb_fasta_seqs_file = fasta_seqs_dir / f'camp_pos_pdb_fasta_seqs.csv'
gold_human_mid_len_pos_pdb_fasta_seqs_natural_aa_file = fasta_seqs_dir / 'positive112_natural_aa_human_200_300mer.csv'
camp_pos_pdb_fasta_seqs_uniprot_file = fasta_seqs_dir / f'camp_pos_pdb_fasta_seqs_uniprot.csv'


def create_camp_positive_pdb_seqs(id=0, check_natural=True):
    """
    Returns:
        fasta file
    """
    unique_positive_data = pd.read_csv(pos_only_seqs_pair_df_file)
    ic(unique_positive_data.columns)
    pdb_id, pep_chain, prot_chain, pep_seq, prot_seq = unique_positive_data.iloc[id].tolist()

    out_file = out_dir / f'camp_{pdb_id}_{prot_chain}_{pep_chain}_pos.fasta'
    save_seq_pair(check_natural, prot_chain, pep_chain, pdb_id, pdb_id, out_file, pep_seq, prot_seq)


def create_camp_positive_pdb_fasta_seqs(id=0, check_natural=True):
    """
    Returns:
        fasta file
    """
    unique_positive_data = pd.read_csv(camp_pos_pdb_fasta_seqs_file)
    ic(unique_positive_data.columns)
    row = unique_positive_data.iloc[id]
    pdb_id, pep_chain, prot_chain, pep_seq, prot_seq, pep_pdb_fasta_seq, prot_pdb_fasta_seq = row.tolist()

    out_file = out_dir / f'camp_{pdb_id}_{prot_chain}_{pep_chain}_pdb_fasta_pos.fasta'
    save_seq_pair(check_natural, prot_chain, pep_chain, pdb_id, pdb_id, out_file, pep_pdb_fasta_seq, prot_pdb_fasta_seq)


def create_camp_positive_uniprot_protein_seqs(id=0, check_natural=True):
    """  """
    unique_positive_data = pd.read_csv(camp_pos_pdb_fasta_seqs_uniprot_file)
    ic(unique_positive_data.columns)
    row = unique_positive_data.iloc[id]
    pdb_id, pep_chain, prot_chain, pep_seq, prot_seq, pep_pdb_fasta_seq, prot_pdb_fasta_seq, uniprot_seq = row.tolist()

    out_file = out_dir / f'camp_{pdb_id}_{prot_chain}_{pep_chain}_pdb_fasta_uniprot_pos.fasta'
    save_seq_pair(check_natural, prot_chain, pep_chain, pdb_id, pdb_id, out_file, pep_pdb_fasta_seq, uniprot_seq)


def save_seq_pair(check_natural, prot_chain, pep_chain, pdb_id_prot, pdb_id_pep, out_file:Path, pep_seq, prot_seq):
    """ prot seq is first """
    if out_file.is_file():
        logger.info(f'out_file {out_file} exits, not overwite')
        return
    if check_natural:
        is_peq_seq_natural, is_prot_seq_natural = check_seq_pair_natural(pep_seq, prot_seq)
    else:
        is_peq_seq_natural, is_prot_seq_natural = True, True
    if is_peq_seq_natural and is_prot_seq_natural:
        logger.info('all natural aa in both seqs')
        ic(out_file.name)
        with open(out_file, 'w', encoding='utf-8') as f:
            f.write(f'>{pdb_id_prot}_{prot_chain}\n{prot_seq}\n')
            f.write(f'>{pdb_id_pep}_{pep_chain}\n{pep_seq}\n')


def check_seq_pair_natural(pep_seq, prot_seq):
    not_natural_aas_pep = get_not_natural_aas_only_supper(pep_seq)
    not_natural_aas_prot = get_not_natural_aas_only_supper(prot_seq)
    is_peq_seq_natural = len(not_natural_aas_pep) == 0
    is_prot_seq_natural = len(not_natural_aas_prot) == 0
    logger.info(f'is_peq_seq_natural {is_peq_seq_natural}')
    logger.info(f'is_prot_seq_natural {is_prot_seq_natural}')
    if not is_peq_seq_natural:
        logger.info(f'is_peq_seq_natural not\n{not_natural_aas_pep}\n{pep_seq}')
    if not is_prot_seq_natural:
        logger.info(f'is_prot_seq_natural not\n{not_natural_aas_prot}\n{prot_seq}')
    return is_peq_seq_natural, is_prot_seq_natural


def create_gold_human_pdb_fasta_seq_input(number=2, select_positive=True, check_natural=True):
    """  Actually only postive, column names: prot_pdb_fasta_seq, pep_pdb_fasta_seq.
    
    Seq: prot first.
    """
    file = gold_human_mid_len_pos_pdb_fasta_seqs_natural_aa_file

    # pos, 'pdb_id', 'prot_chain', 'pep_chain'
    df = pd.read_csv(file)
    ic(df.head(1))
    ic(df.columns)
    ic(len(df))

    for id in range(number):
        row = df.iloc[id]
        prot_chain = row['prot_chain']
        pep_chain = row['pep_chain']
        if select_positive:
            pdb_id_prot = row['pdb_id']
            pdb_id_pep = pdb_id_prot
            out_file = out_dir / f'{pdb_id_prot}_{prot_chain}_{pep_chain}_pdb_fasta_pos_human.fasta'
        else:
            pdb_id_prot = row['pdb_id_prot']
            pdb_id_pep = row['pdb_id_pep']
            out_file = out_dir / f'{pdb_id_prot}_{prot_chain}_{pdb_id_pep}_{pep_chain}_pdb_fasta_neg_human.fasta'
        pep_seq = row['pep_pdb_fasta_seq']
        prot_seq = row['prot_pdb_fasta_seq']
        logger.info(f'{pdb_id_prot}_{prot_chain}_{pdb_id_pep}_{pep_chain}')
        save_seq_pair(check_natural, prot_chain, pep_chain, pdb_id_prot, pdb_id_pep, out_file, pep_seq, prot_seq)


def create_gold_pdb_fasta_seq_input(id=0, select_positive=True, check_natural=True):
    """  column names: prot_pdb_fasta_seq, pep_pdb_fasta_seq.

    Seq: prot first.
    """
    if select_positive:
        file = gold_pos_pdb_fasta_seqs_natural_aa_file
    else:
        file = gold_neg_pdb_fasta_seqs_natural_aa_file

    # pos, 'pdb_id', 'prot_chain', 'pep_chain'
    # neg, 'pdb_id_prot', 'prot_chain', 'pdb_id_pep', 'pep_chain'
    df = pd.read_csv(file)
    ic(df.head())

    row = df.iloc[id]
    prot_chain = row['prot_chain']
    pep_chain = row['pep_chain']
    if select_positive:
        pdb_id_prot = row['pdb_id']
        pdb_id_pep = pdb_id_prot
        out_file = out_dir / f'{pdb_id_prot}_{prot_chain}_{pep_chain}_pdb_fasta_pos.fasta'
    else:
        pdb_id_prot = row['pdb_id_prot']
        pdb_id_pep = row['pdb_id_pep']
        out_file = out_dir / f'{pdb_id_prot}_{prot_chain}_{pdb_id_pep}_{pep_chain}_pdb_fasta_neg.fasta'
    pep_seq = row['pep_pdb_fasta_seq']
    prot_seq = row['prot_pdb_fasta_seq']

    save_seq_pair(check_natural, prot_chain, pep_chain, pdb_id_prot, pdb_id_pep, out_file, pep_seq, prot_seq)


def create_gold_pdb_seq_input(pdb_id_chain: str, positive=None, check_natural=True):
    """  """
    items = pdb_id_chain.split('_')
    if len(items) == 3:
        pdb_id_prot, prot_chain, pep_chain = items
        pdb_id_pep = pdb_id_prot
        if positive is None:
            positive = True
    elif len(items) == 4:
        pdb_id_prot, prot_chain, pdb_id_pep, pep_chain = items
        if positive is None:
            positive = False
    prot_seq = get_pdb_seq_from_pdb_file(pdb_id_prot, prot_chain)
    pep_seq = get_pdb_seq_from_pdb_file(pdb_id_pep, pep_chain)
    if positive:
        if len(items) == 3:
            out_file = out_dir / f'{pdb_id_prot}_{prot_chain}_{pep_chain}_pos.fasta'
        else:
            out_file = out_dir / f'{pdb_id_prot}_{prot_chain}_{pdb_id_pep}_{pep_chain}_pos.fasta'
    else:
        out_file = out_dir / f'{pdb_id_prot}_{prot_chain}_{pdb_id_pep}_{pep_chain}_neg.fasta'
    save_seq_pair(check_natural, prot_chain, pep_chain, pdb_id_prot, pdb_id_pep, out_file, pep_seq, prot_seq)


def batch_create_gold_pdb_seq_input():
    """  """
    pdb_id_chains = [
        # '6peu_A_M',
        # '6s6q_A_6mf6_C',
        # '6i51_H_I',
        # '6itm_A_B',
        '5qtu_A_6kac_V',
    ]
    for pdb_id_chain in pdb_id_chains:
        create_gold_pdb_seq_input(pdb_id_chain)


def create_create_camp_positive(id=2):
    """  """
    create_camp_positive_pdb_seqs(id)
    create_camp_positive_pdb_fasta_seqs(id)
    create_camp_positive_uniprot_protein_seqs(id)

if __name__ == "__main__":
    # create_gold_pdb_fasta_seq_input(id=0, select_positive=True)
    # create_gold_pdb_fasta_seq_input(id=1, select_positive=0, check_natural=True)
    batch_create_gold_pdb_seq_input()

    # create_gold_human_pdb_fasta_seq_input()
    
    pass