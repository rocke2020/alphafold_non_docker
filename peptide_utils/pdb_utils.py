from pathlib import Path
import pandas as pd
from pandas import DataFrame
import os, sys, re, json
sys.path.append(os.path.abspath('.'))
from utils.log_util import logger
from tqdm import tqdm
from peptide_utils.pep_utils import is_natural_only_supper, aminoacids
from peptide_utils.proteins.shared.protein import get_seq_from_pdb_file_pdb_parser
from peptide_utils.aa_utils import basic_aa_3chars_to_1chars
from typing import List


pdb_files_dir = Path('data/pdb/pdb_files')
pdb_fasta_seqres_file = Path('data/pdb/pdb_seqres.txt')
pdb_fasta_seqres_dict_file = Path('data/pdb/pdb_seqres.json')
organism_prefixes = ('ORGANISM_SCIENTIFIC', 'ORGANISM_COMMON')
# TODO maybe only Homo sapiens, because Pediculus humanus is 头虱/人虱.
human_organism_names = ('HOMO SAPIENS', 'HUMAN')
HUMAN_ORGANISM = 'homo sapiens'
PDB_ID = 'pdb_id'
PEP_CHAIN = 'pep_chain'
PROT_CHAIN = 'prot_chain'
PEP_SEQ = 'pep_seq'
PROT_SEQ = 'prot_seq'
UNIPROT_ID = 'Uniprot_id'
PROTEIN_FAMILIES = 'protein_families'
SEQUENCE = 'Sequence'

# in real usage read the pdb_fasta_seqres_dict_file
# with open(pdb_fasta_seqres_dict_file, 'r', encoding='utf-8') as f:
    # pdb_fasta_seqres_dict = json.load(f)
pdb_fasta_seqres_dict = {}


def natural_aa_ratio(seq):
    """  """
    total_len = len(seq)
    count = 0
    for aa in seq:
        if aa in aminoacids:
            count += 1
    return count / total_len


def get_pdb_seq_from_pdb_file(pdb_id, seq_chain, pdb_files_dir=pdb_files_dir):
    """
    Special case:
        if seq_chain == 'na': return 'NA'
    """
    pdb_files_dir = Path(pdb_files_dir)
    if seq_chain == 'na':
        return 'NA'
    pdb_file = pdb_files_dir / f'{pdb_id}.pdb'
    return get_seq_from_pdb_file_pdb_parser(pdb_file, seq_chain)


def get_pdb_fasta_seq_from_dict(pdb_id, chain):
    """  """
    pdb_id_chain = pdb_id.lower() + '_' + chain
    seq  = pdb_fasta_seqres_dict.get(pdb_id_chain, '')
    if seq == '':
        logger.info(f'pdb_id_chain {pdb_id_chain} has not fasta seq in dict')
    return seq


def write_pdb_fasta_seqs_to_dict():
    """  """
    pdb_id_chain_seqs = {}
    pdb_id_line_head = re.compile(r'^>\d[\da-zA-Z]{3}_[\da-zA-Z]+')
    seq_count = 0
    with open(pdb_fasta_seqres_file,'r') as f:
        for line in f.readlines():
            line = line.strip()
            if pdb_id_line_head.search(line):
                pdb_id = line[1:5]
                items = line.split()
                assert items[1].startswith('mol:')
                chain = items[0][6:]
                seq_count = 0
                pdb_id_chain = pdb_id + '_' + chain
            else:
                seq_count += 1
                if seq_count > 1:
                    logger.error(f'Error line for continous {line}')
                    return
                pdb_id_chain_seqs[pdb_id_chain] = line
    with open(pdb_fasta_seqres_dict_file, 'w', encoding='utf-8') as f:
        json.dump(pdb_id_chain_seqs, f, ensure_ascii=False, indent=4)


def read_pepbdb_seq(file:Path, seq_chain=None):
    """ specialy dataset pepbdb with binding sites """
    aas = []
    with open(file, 'r', encoding='utf-8') as f:
        last_index = ''
        for line in f:
            items = line.split()
            if items[0] == 'ATOM':
                amino_acid = items[3][-3:]
                aa = basic_aa_3chars_to_1chars.get(amino_acid, 'X')
                chain = items[4]
                if seq_chain:
                    pdb_id = file.parent.name.split('_')[0]
                    index = get_index_with_seq_chain(pdb_id, seq_chain, items, chain)
                else:
                    if len(chain) == 1:
                        index = items[5]
                    elif len(chain) > 1:
                        index = chain[1:]
                    else:
                        logger.info(f'{file.parent.name} has abnormal chain at line {items}')
                        raise Exception()
                if index != last_index:
                    aas.append(aa)
                last_index = index
            elif items[0] == 'TER':
                break
            elif items[0] == 'HETATM':
                logger.debug('items[0] is abnormal: HETATM, and so skip it')
                continue
    return ''.join(aas)


def get_index_with_seq_chain(pdb_id, seq_chain, items, chain):
    if chain == seq_chain:
        index = items[5]
    elif chain[0] == seq_chain:
        index = chain[1:]
    else:
        logger.info(f'{pdb_id} has abnormal chain at line {seq_chain}, items: {items}')
        raise Exception()
    return index


def get_all_downloaded_pdb_files(pdb_files_dir=pdb_files_dir):
    """  """
    all_downloaded_pdb_ids = []
    for file in pdb_files_dir.glob('*.pdb'):
        all_downloaded_pdb_ids.append(file.stem.lower())
    return all_downloaded_pdb_ids


def write_undownloaded_ids_to_file(un_downloaded_pdb_ids, out_file='data_process/batch_download_pdb/list_file.txt'):
    """  """
    out_str = ','.join(un_downloaded_pdb_ids) + '\n'
    with open(out_file, 'w', encoding='utf-8') as f:
        f.write(out_str)


def create_seqs_from_pdb_id_and_chain(
    id_and_chains, only_natural_in_protein=True, check_id='6qmp', peptide_length_filter=True):
    """  """
    items = []
    for pdb_id, pep_chain, prot_chain in tqdm(id_and_chains):
        try:
            pdb_pep_seq = get_pdb_seq_from_pdb_file(pdb_id, pep_chain)
            pdb_prot_seq = get_pdb_seq_from_pdb_file(pdb_id, prot_chain)
            items.append([pdb_id, pep_chain, prot_chain, pdb_pep_seq, pdb_prot_seq])
            if check_id and pdb_id == check_id:
                logger.info(f'pdb_id {check_id} {pep_chain} {pdb_pep_seq}')
                logger.info(f'pdb_id {check_id} {prot_chain} {pdb_prot_seq}')
        except Exception as identifier:
            logger.exception(f'{pdb_id} {pep_chain} {prot_chain} error: {identifier}')

    new_df = DataFrame(
        items, columns=[PDB_ID, PEP_CHAIN, PROT_CHAIN, PEP_SEQ, PROT_SEQ])
    normal_peptide_aa_df = new_df[new_df[PEP_SEQ].map(is_natural_only_supper)]
    logger.info(f'peptide_normal_aa_df.shape {normal_peptide_aa_df.shape}')
    if only_natural_in_protein:
        seq_pairs_df = normal_peptide_aa_df[normal_peptide_aa_df[PROT_SEQ].map(is_natural_only_supper)]
        logger.info(f'after protein normal_aa, df.shape {seq_pairs_df.shape}')
    # There are duplicate seq pairs while the pdb-id-chains are different.
    seq_pairs_df = seq_pairs_df.drop_duplicates([PEP_SEQ, PROT_SEQ])
    seq_pairs_df = seq_pairs_df.dropna(subset=[PEP_SEQ, PROT_SEQ])
    logger.info(f'after drop_duplicates seqs_pair_df shape: {seq_pairs_df.shape}')
    seq_pairs_df = seq_pairs_df[seq_pairs_df[PROT_SEQ].map(lambda x: 50<len(x)<5000)]
    if peptide_length_filter:
        seq_pairs_df = seq_pairs_df[seq_pairs_df[PEP_SEQ].map(lambda x: 5<=len(x)<=50)]
        logger.info(f'after peptide and protein lenght filter, shape: {seq_pairs_df.shape}')
    return seq_pairs_df


def is_hunman_type_by_checking_pdb(pdb_id, pdb_files_dir=pdb_files_dir):
    """  """
    pdb_file = pdb_files_dir / f'{pdb_id}.pdb'
    organisms = []
    with open(pdb_file, 'r', encoding='utf-8') as f:
        for line in f:
            for organism_prefix in organism_prefixes:
                if organism_prefix in line:
                    items = line.strip().split(':')
                    organism_type = items[-1]
                    organisms.append(organism_type)
                    break
    for organism in organisms:
        for human_organism_name in human_organism_names:
            if human_organism_name in organism:
                return True
    return False


def move_pdb_file_to_tmp_folder_for_copy(pdb_id):
    """  """
    pass


if __name__ == "__main__":
    write_pdb_fasta_seqs_to_dict()
    pass