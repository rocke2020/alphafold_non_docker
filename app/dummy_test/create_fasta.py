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


def get_the_protein_pdb_seq(id=0):
    """  
    Returns:
        fasta file
    """
    unique_positive_data = pd.read_csv(pos_only_seqs_pair_df_file)
    ic(unique_positive_data.columns)
    pdb_id, pep_chain, prot_chain, pep_seq, prot_seq = unique_positive_data.iloc[id].tolist()
    out_file = out_dir / f'{pdb_id}_{prot_chain}.fasta'
    with open(out_file, 'w', encoding='utf-8') as f:
        f.write(f'>{pdb_id}_{prot_chain}\n{prot_seq}\n')

get_the_protein_pdb_seq()