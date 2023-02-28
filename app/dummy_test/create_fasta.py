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
gold_test_file = test_data_dir / 'test-complex-num112-all positive.CSV'


def get_the_protein_pdb_seq_camp_all_positive(id=1):
    """  
    Returns:
        fasta file
    """
    unique_positive_data = pd.read_csv(pos_only_seqs_pair_df_file)
    ic(unique_positive_data.columns)
    pdb_id, pep_chain, prot_chain, pep_seq, prot_seq = unique_positive_data.iloc[id].tolist()
    out_file = out_dir / f'{pdb_id}_{pep_chain}_{prot_chain}.fasta'
    ic(out_file.name)
    with open(out_file, 'w', encoding='utf-8') as f:
        f.write(f'>{pdb_id}_{pep_chain}\n{pep_seq}\n')
        f.write(f'>{pdb_id}_{prot_chain}\n{prot_seq}\n')


get_the_protein_pdb_seq_camp_all_positive()