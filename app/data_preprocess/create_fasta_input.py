import re, random
from pathlib import Path
import json
import pandas as pd
import numpy as np
from pandas import DataFrame
import os, sys
from icecream import ic
ic.configureOutput(includeContext=True, argToStringFunction=lambda _: str(_))


version='2.3'
leiyipin_dir = Path('/mnt/sda/bio_drug_corpus/ppi/ppi-leiyipin/peptide_protein_interaction_v2')
pos_only_seqs_pair_df_file = leiyipin_dir / f'peptide_protein_interaction_v{version}_unique.csv'

out_dir = Path('/mnt/sdc/af_input')
test_data_dir = Path('/mnt/sdb/bio_drug_corpus/ppi/Improving Peptide-Protein Docking with AlphaFold-Multimer using Forced Sampling')
gold_positive_test_file = test_data_dir / 'test-complex-num112-all positive.CSV'
gold_negative_test_file = test_data_dir / 'test-complex-num672-1pos5neg.CSV'

test_data_dir = Path('/mnt/sdb/bio_drug_corpus/ppi/Improving Peptide-Protein Docking with AlphaFold-Multimer using Forced Sampling')
pdb_seqs_dir = test_data_dir / 'pdb_seqs'
gold_pos_test_pdb_seqs_file = pdb_seqs_dir / 'positive112.csv'
gold_neg_test_pdb_seqs_file = pdb_seqs_dir / 'negative560.csv'
fasta_seqs_dir = test_data_dir / 'fasta_seqs'
gold_pos_test_pdb_fasta_seqs_file = fasta_seqs_dir / 'positive112.csv'
gold_neg_test_pdb_fasta_seqs_file = fasta_seqs_dir / 'negative560.csv'


def create_camp_all_positive_input(id=1):
    """
    Returns:
        fasta file
    """
    unique_positive_data = pd.read_csv(pos_only_seqs_pair_df_file)
    ic(unique_positive_data.columns)
    pdb_id, pep_chain, prot_chain, pep_seq, prot_seq = unique_positive_data.iloc[id].tolist()
    out_file = out_dir / f'{pdb_id}_{prot_chain}_{pep_chain}.fasta'
    ic(out_file.name)
    with open(out_file, 'w', encoding='utf-8') as f:
        f.write(f'>{pdb_id}_{prot_chain}\n{prot_seq}\n')
        f.write(f'>{pdb_id}_{pep_chain}\n{pep_seq}\n')


def create_gold_pdb_fasta_seq_input(id=0, select_positive=True):
    """  column names: prot_pdb_fasta_seq, pep_pdb_fasta_seq.

    Seq: prot first.
    """
    if select_positive:
        file = gold_pos_test_pdb_fasta_seqs_file
    else:
        file = gold_neg_test_pdb_fasta_seqs_file

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
    ic(out_file.name)
    with open(out_file, 'w', encoding='utf-8') as f:
        f.write(f'>{pdb_id_prot}_{prot_chain}\n{prot_seq}\n')
        f.write(f'>{pdb_id_pep}_{pep_chain}\n{pep_seq}\n')


if __name__ == "__main__":
    create_gold_pdb_fasta_seq_input(id=0, select_positive=True)
    create_gold_pdb_fasta_seq_input(id=0, select_positive=0)
    create_camp_all_positive_input()