import re, random
from pathlib import Path
import pickle
import json
import pandas as pd
import numpy as np
from pandas import Series, DataFrame
import os, sys, shutil
from icecream import ic
ic.configureOutput(includeContext=True, argToStringFunction=lambda _: str(_))
sys.path.append(os.path.abspath('.'))
from utils.log_util import logger
from functools import partial
from collections import defaultdict
from peptide_utils.pdb_utils import (
    PDB_ID, PEP_CHAIN, PEP_SEQ, PROT_SEQ, SEQUENCE, UNIPROT_ID, HUMAN_ORGANISM, PROTEIN_FAMILIES,
)


pdb_root_dir = Path('data/pdb')
uniport_root_dir = Path('data/uniprot')
raw_uniprot2seq_file = uniport_root_dir / 'raw_uniprot2seq.csv'
raw_uniprot2OS_file = uniport_root_dir / 'raw_uniprot2OS.csv'
uniprot2seq_and_prot_families_file = uniport_root_dir / 'uniprot2seq.tsv'
uniport_id_list_file = 'data/uniport/uniport_id_list.tsv'
pdb_chain_uniprot_file = pdb_root_dir / 'pdb_chain_uniprot.tsv'
unique_pdb_chain_uniprot_file = pdb_root_dir / 'unique_pdb_chain_uniprot.csv'
pdb_id_and_chain_to_uniprot_id_dict_file = pdb_root_dir / 'pdb_id_and_chain_to_uniprot_id_dict.json'
uniprot_sprot_fasta_file = uniport_root_dir / 'uniprot_sprot.fasta'

similar_file = 'data/uniprot/similar'
uniprot2protein_families_file = uniport_root_dir / 'uniprot2protein_families.csv'  # out file

PDB = 'PDB'
CHAIN = 'CHAIN'
SP_PRIMARY = 'SP_PRIMARY'
SP_BEG = 'SP_BEG'
SP_END = 'SP_END'
UNIPROT_SP_FRAGMENT_SEQ = 'uniprot_sp_fragment_seq'


def read_uniprot_sprot_fasta(read_cache=True):
    """
    uniprot_sprot_fasta fragments:
        >sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-001R PE=4 SV=1
        MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPS
        EKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLD
        AKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHL
        EKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDD
        SFRKIYTDLGWKFTPL
        >sp|Q6GZX3|002L_FRG3G Uncharacterized protein 002L OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-002L PE=4 SV=1
        MSIIGATRLQNDKSDTYSAGPCYAGGCSAFTPRGTCGKDWDLGEQTCASGFCTSQPLCAR
    """
    check_id = 'Q6GZX4'
    if read_cache and raw_uniprot2seq_file.is_file():
        with open(raw_uniprot2seq_file, 'rb') as f:
            uniprot_sp_id2seq = pickle.load(f)
        with open(raw_uniprot2OS_file, 'rb') as f:
            uniprot_sp_id2OS = pickle.load(f)
        # for check
        seq = uniprot_sp_id2seq[check_id]
        logger.info(f'sp id: {check_id}, seq {seq}')
        _os = uniprot_sp_id2OS[check_id]
        logger.info(f'len(_os) {len(_os)}, os {_os}')
        return uniprot_sp_id2seq, uniprot_sp_id2OS

    raw_uniprot_sp_id2seq = defaultdict(list)
    uniprot_sp_id2OS = dict()
    with open(uniprot_sprot_fasta_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>sp|'):
                items = line.split('|')
                sp_primary = items[1]
                post_info = items[2]
                start_i = post_info.index('OS=') + 3
                end_i = post_info.index('OX=')-1
                sp_os = post_info[start_i: end_i].lower()
                uniprot_sp_id2OS[sp_primary] = sp_os
            else:
                raw_uniprot_sp_id2seq[sp_primary].append(line)

    uniprot_sp_id2seq = {k: ''.join(v) for k, v in raw_uniprot_sp_id2seq.items()}
    ic(len(uniprot_sp_id2seq))

    # for check
    seq = uniprot_sp_id2seq[check_id]
    logger.info(f'sp id: {check_id}, seq {seq}')
    _os = uniprot_sp_id2OS[check_id]
    logger.info(f'len(_os) {len(_os)}, os {_os}')
    with open(raw_uniprot2seq_file, 'wb') as f:
        pickle.dump(uniprot_sp_id2seq, f)
    with open(raw_uniprot2OS_file, 'wb') as f:
        pickle.dump(uniprot_sp_id2OS, f)
    return uniprot_sp_id2seq, uniprot_sp_id2OS


uniprot_sp_id2seq, uniprot_sp_id2OS = read_uniprot_sprot_fasta()
with open(pdb_id_and_chain_to_uniprot_id_dict_file, 'r', encoding='utf-8') as f:
    pdb_id_and_chain_to_uniprot_id_dict = json.load(f)


def get_sp_range(row):
    """  """
    beg = row[SP_BEG]
    end = row[SP_END]
    sp_range = f'{beg}-{end}'
    return sp_range


def anaylyze_pdb_chain_uniprot():
    """
    https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
    df_sifts from pdb_chain_uniprot_file:
        # 2022/08/08 - 11:12 | PDB: 31.22 | UniProt: 2022.03
        PDB	CHAIN	SP_PRIMARY	RES_BEG	RES_END	PDB_BEG	PDB_END	SP_BEG	SP_END
        101m	A	P02185	1	154	0	153	1	154
        102l	A	P00720	1	40	1	40	1	40
        102l	A	P00720	42	165	41	None	41	164
        102m	A	P02185	1	154	0	153	1	154
    sum:
        1. in the pdb_chain_uniprot, there are duplicate SP_PRIMARY
    """
    pdb_chain_uniprot_file = 'data/pdb/pdb_chain_uniprot.tsv'
    df = pd.read_csv(pdb_chain_uniprot_file, skiprows=1, sep='\t')
    ic(df.columns.tolist())

    df['count'] = df.apply(get_sp_range, axis=1)
    ic(df['count'].head())

    _df = df.groupby(SP_PRIMARY).agg({'count': lambda x: len(set(x))})
    _df = _df.reset_index()
    ic(_df.head())
    _df = _df.sort_values(by=['count'], ascending=False)
    ic(_df.head())
    count = _df['count'].max()
    ic(count)


def choose_item_with_SP_Primiary_index(row, column_name):
    """ TODO to replace """
    index = -1
    for i, item in enumerate(row[SP_PRIMARY]):
        if item.startswith('P'):
            index = i
            break
    item = list(row[column_name])[index]
    return item


def choose_SP_Primiary_starts_with_P(array):
    """ TODO it is improper to select choose_SP_Primiary_starts_with_P, but prefer to select the longest sp full seq.
    """
    for item in array:
        if item.startswith('P'):
            return item
    return item


def choose_SP_Primiary(array):
    """  """
    proper_item = ''
    max_len = -1
    for item in array:
        seq_os = uniprot_sp_id2OS.get(item, '')
        if HUMAN_ORGANISM in seq_os:
            return item
        seq_len = len(uniprot_sp_id2seq.get(item, ''))
        if seq_len > max_len:
            max_len = seq_len
            proper_item = item
    return proper_item


def parse_pdb_chain_uniprot(read_cache=True):
    """ let pdb
    https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
    df_sifts from pdb_chain_uniprot_file:
        # 2022/08/08 - 11:12 | PDB: 31.22 | UniProt: 2022.03
        PDB	CHAIN	SP_PRIMARY	RES_BEG	RES_END	PDB_BEG	PDB_END	SP_BEG	SP_END
        101m	A	P02185	1	154	0	153	1	154
        102l	A	P00720	1	40	1	40	1	40
        102l	A	P00720	42	165	41	None	41	164
        102m	A	P02185	1	154	0	153	1	154

    one pdb-id and chain links to different uniProt Primary id
        duplicated_df_sifts.head():
        PDB	CHAIN	SP_PRIMARY	RES_BEG	RES_END	PDB_BEG	PDB_END	SP_BEG	SP_END
        1ao7	D	A0A075B6T6	1	90	1	92	23	112
        1ao7	D	P01848	112	204	118	None	1	93
        8dfl	D	P22001	1	575	None	None	1	575
        8dfl	D	P42212	590	826	None	None	2	238
    
    Returns:
        df, columns: [PDB, CHAIN, SP_PRIMARY]
    """
    if read_cache and unique_pdb_chain_uniprot_file.is_file():
        df_sifts = pd.read_csv(unique_pdb_chain_uniprot_file)
        ic(len(df_sifts))
        if len(df_sifts) > 0:
            return df_sifts
    raw_df_sifts = pd.read_csv(pdb_chain_uniprot_file, sep = '\t', skiprows=1)
    raw_df_sifts = raw_df_sifts[[PDB, CHAIN, SP_PRIMARY]]
    logger.info(f'raw_df_sifts len {len(raw_df_sifts)}')

    # Use df_sifts_del to pre view what will be filtered by df_sifts[CHAIN].isna
    df_sifts_chain_isna = raw_df_sifts[raw_df_sifts[CHAIN].isna()]
    logger.info(f'df_sifts chain isna len {len(df_sifts_chain_isna)}')
    logger.info(f'df_sifts chain isna\n{df_sifts_chain_isna.head()}')

    full_df_sifts = raw_df_sifts[~raw_df_sifts[CHAIN].isna()]
    logger.info(f'full_df_sifts len {len(full_df_sifts)}')

    df_sifts = full_df_sifts.drop_duplicates([PDB, CHAIN, SP_PRIMARY])
    logger.info(f'df_sifts len {len(df_sifts)}')

    # When 1 pdb-id and chain links to multiple uniProt Primary, choose the one uniprot id
    duplicated_df_sifts = df_sifts[df_sifts.duplicated([PDB, CHAIN], keep=False)]
    ic(len(duplicated_df_sifts))
    ic(duplicated_df_sifts.head())

    duplicated_df_sifts = duplicated_df_sifts.groupby(
        by=[PDB, CHAIN]).agg({SP_PRIMARY: choose_SP_Primiary})
    ic(len(duplicated_df_sifts))
    duplicated_df_sifts = duplicated_df_sifts.reset_index()
    ic(duplicated_df_sifts.head())

    _df_sifts = df_sifts.drop_duplicates([PDB, CHAIN], keep=False)
    ic(len(_df_sifts))
    df_sifts = pd.concat([_df_sifts, duplicated_df_sifts], axis=0)
    ic(len(df_sifts))  # 622914
    df_sifts.to_csv(unique_pdb_chain_uniprot_file, index=False, sep=',')
    return df_sifts


def parse_pdb_chain_uniprot_seq_fragments():
    """ TODO not finished but choose to select full uniprot seq.
    parse_pdb_chain_uniprot_with_seq is much more complex than parse_pdb_chain_uniprot(without seq) and so make 2 funcs.
    # to select the longest sp full seq.
    https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
    df_sifts from pdb_chain_uniprot_file:
        # 2022/08/08 - 11:12 | PDB: 31.22 | UniProt: 2022.03
        PDB	CHAIN	SP_PRIMARY	RES_BEG	RES_END	PDB_BEG	PDB_END	SP_BEG	SP_END
        101m	A	P02185	1	154	0	153	1	154
        102l	A	P00720	1	40	1	40	1	40
        102l	A	P00720	42	165	41	None	41	164
        102m	A	P02185	1	154	0	153	1	154

    one pdb-id and chain links to different uniProt Primary id
        duplicated_df_sifts.head():
        PDB	CHAIN	SP_PRIMARY	RES_BEG	RES_END	PDB_BEG	PDB_END	SP_BEG	SP_END
        1ao7	D	A0A075B6T6	1	90	1	92	23	112
        1ao7	D	P01848	112	204	118	None	1	93
        8dfl	D	P22001	1	575	None	None	1	575
        8dfl	D	P42212	590	826	None	None	2	238
    """
    raw_df_sifts = pd.read_csv(pdb_chain_uniprot_file, sep = '\t', skiprows=1)
    raw_df_sifts = raw_df_sifts[[PDB, CHAIN, SP_PRIMARY, SP_BEG, SP_END]]
    logger.info(f'raw_df_sifts len {len(raw_df_sifts)}')

    # Use df_sifts_del to pre view what will be filtered by df_sifts[CHAIN].isna
    df_sifts_del = raw_df_sifts[raw_df_sifts[CHAIN].isna()]
    logger.info(f'df_sifts_del len {len(df_sifts_del)}')
    logger.info(f'df_sifts_del\n{df_sifts_del.head()}')

    full_df_sifts = raw_df_sifts[~raw_df_sifts[CHAIN].isna()]
    logger.info(f'full_df_sifts len {len(full_df_sifts)}')

    # with_seq True, firstly del the three dulicates [PDB, CHAIN, SP_PRIMARY], then filters only [PDB, CHAIN]
    # duplicate, and finally adds back the three dulicates.
    df_single_uniprot_fragment = full_df_sifts.drop_duplicates([PDB, CHAIN, SP_PRIMARY], keep=False)
    logger.info(f'df_sifts len {len(df_single_uniprot_fragment)}')

    # When 1 pdb-id and chain links to multiple uniProt Primary, choose the one starts with 'P' if 'P' is available.
    duplicated_df_sifts = df_single_uniprot_fragment[df_single_uniprot_fragment.duplicated([PDB, CHAIN], keep=False)]
    ic(len(duplicated_df_sifts))
    ic(duplicated_df_sifts.head())

    duplicated_df_sifts = duplicated_df_sifts.groupby(by=[PDB, CHAIN]).agg(lambda x: list(x))
    ic(len(duplicated_df_sifts))
    # TODO select the longest sp full seq.
    duplicated_df_sifts['_SP_PRIMARY'] = duplicated_df_sifts[SP_PRIMARY].map(choose_SP_Primiary_starts_with_P)
    choose_func_SP_BEG = partial(choose_item_with_SP_Primiary_index, column_name=SP_BEG)
    choose_func_SP_END = partial(choose_item_with_SP_Primiary_index, column_name=SP_END)
    duplicated_df_sifts['_SP_BEG'] = duplicated_df_sifts.apply(choose_func_SP_BEG, axis=1)
    duplicated_df_sifts['_SP_END'] = duplicated_df_sifts.apply(choose_func_SP_END, axis=1)
    duplicated_df_sifts.drop(columns=[SP_PRIMARY, SP_BEG, SP_END], inplace=True)
    duplicated_df_sifts.rename(columns={'_SP_PRIMARY': SP_PRIMARY, '_SP_BEG': SP_BEG, '_SP_END':SP_END}, inplace=True)
    duplicated_df_sifts = duplicated_df_sifts.reset_index()
    ic(duplicated_df_sifts.head())

    _df_sifts = df_single_uniprot_fragment.drop_duplicates([PDB, CHAIN], keep=False)
    ic(len(_df_sifts))
    df_single_uniprot_fragment = pd.concat([_df_sifts, duplicated_df_sifts], axis=0)
    ic(len(df_single_uniprot_fragment))

    df = add_uniprot_seq_fragments_to_df(df_single_uniprot_fragment, full_df_sifts)
    ic(len(df))  # 622914
    return df


def add_uniprot_seq_fragments_to_df(df_single_uniprot_fragment, full_df_sifts):
    """
    fragments_df is pdb-chain links to multiple SP_PRIMARY sequence fragments
        PDB	CHAIN	SP_PRIMARY	RES_BEG	RES_END	PDB_BEG	PDB_END	SP_BEG	SP_END
       102l	A	P00720	1	40	1	40	1	40
       102l	A	P00720	42	165	41	None	41	164
       102l	A	P00721	1	40	1	40	1	40
       102l	A	P00721	42	165	41	None	41	164
        1es0	B	Q05329	7	22	207P	None	207	230
        1es0	B	Q31135	23	24	None	None	24	25
        1es0	B	Q05329	25	28	None	None	233	236
        1es0	B	Q31135	29	213	None	188	30	214

    TODO select only one SP_PRIMARY for mulitp SP_PRIMARY
    When 1 id-chain links 2 SP_PRIMARY
    1. each SP_PRIMARY have multiple fragments
    2. 1 SP_PRIMARY have multiple fragments and the other has 1 fragment, after
    pd.concat([df_single_uniprot_fragment, fragments_df])

    """

    fragments_df = full_df_sifts[full_df_sifts.duplicated([PDB, CHAIN, SP_PRIMARY], keep=False)]
    df = pd.concat([df_single_uniprot_fragment, fragments_df])
    df = df.groupby([PDB, CHAIN]).agg(lambda x: list(x)).reset_index()
    ic(len(df))
    ic(df.head())

    uniprot_sp_seqs = []
    for i, row in df.iterrows():
        begin_indices = row[SP_BEG]
        end_indices = row[SP_END]
        ids = row[SP_PRIMARY]
        if i < 5:
            logger.info(f'i {i}, row\n{row}')
        if ids[0] not in uniprot_sp_id2seq:
            seq=''
        else:
            sub_seqs = []
            for beg_i, end_i, sp_id in zip(begin_indices, end_indices, ids):
                try:
                    beg_i = int(beg_i)
                    end_i = int(end_i)
                    sub_seqs.append(uniprot_sp_id2seq[sp_id][beg_i: end_i])
                except Exception as identifier:
                    logger.exception(f'i {i}, row\n{row}\n{identifier}')
                    sub_seqs = []
                    break
            seq = ''.join(sub_seqs)
        uniprot_sp_seqs.append(seq)
    df[UNIPROT_SP_FRAGMENT_SEQ] = uniprot_sp_seqs
    ic(df.head())
    return df


def create_uniprot_to_seq_with_prot_families():
    """
    Returns:
        df, columns: UNIPROT_ID, SEQUENCE, PROTEIN_FAMILIES
    some uniprot id has no sp sequence because it is UniProtKB unreviewed (TrEMBL)
    some id has no protein families
    """
    Uniprot_id = []
    Sequence = []
    uniprot2protein_families, uniprot2protein_families_df = get_protein_families_with_uniport_ids()
    for k, v in uniprot_sp_id2seq.items():
        Uniprot_id.append(k)
        Sequence.append(v)
    seq_df = DataFrame({
        UNIPROT_ID: Uniprot_id,
        SEQUENCE: Sequence,
    })
    ic(len(seq_df))
    ic(seq_df[seq_df[UNIPROT_ID] == 'P00750'])
    uniprot_to_seq_and_prot_families_df = seq_df.merge(uniprot2protein_families_df, how='outer')

    ic(uniprot_to_seq_and_prot_families_df[uniprot_to_seq_and_prot_families_df[UNIPROT_ID] == 'P00750'])
    uniprot_to_seq_and_prot_families_df = uniprot_to_seq_and_prot_families_df.fillna('Unknown_from_uniprot')
    uniprot_to_seq_and_prot_families_df = uniprot_to_seq_and_prot_families_df.drop_duplicates(
        ['Uniprot_id', 'Sequence']).reset_index(drop=True)
    
    ic(uniprot_to_seq_and_prot_families_df.shape)
    ic(uniprot_to_seq_and_prot_families_df.head())
    uniprot_to_seq_and_prot_families_df.to_csv(uniprot2seq_and_prot_families_file, index=False, sep=',')


def get_pdb_id_and_chain_to_uniprot_id_dict():
    """  """
    df_sifts = pd.read_csv(unique_pdb_chain_uniprot_file)
    ic(df_sifts.columns)
    ic(len(df_sifts))
    pdb_id_and_chain_to_uniprot_id_dict = {}
    for i, row in df_sifts.iterrows():
        pdb_id, chain, unprot_id = row.tolist()
        pdb_id_and_chain_to_uniprot_id_dict[f'{pdb_id.lower()}_{chain}'] = unprot_id
    ic(len(pdb_id_and_chain_to_uniprot_id_dict))

    with open(pdb_id_and_chain_to_uniprot_id_dict_file, 'w', encoding='utf-8') as f:
        json.dump(pdb_id_and_chain_to_uniprot_id_dict, f, ensure_ascii=False, indent=4)


def get_uniprot_seq_from_pdb_id_and_chain(pdb_id, chain):
    """  """
    pdb_id_and_chain = f'{pdb_id.lower()}_{chain}'
    uniport_id = pdb_id_and_chain_to_uniprot_id_dict.get(pdb_id_and_chain, None)
    if uniport_id is None:
        logger.info(f'pdb_id_and_chain {pdb_id_and_chain} has not uniprot id')
        return ''
    seq = uniprot_sp_id2seq.get(uniport_id, '')
    if seq == '':
        logger.info(f'pdb_id_and_chain {pdb_id_and_chain} has not uniprot seq')
    return seq


def get_protein_families_with_uniport_ids(save_df=False):
    """  
    Args:
        'GDSL' lipolytic enzyme family. IAH1 subfamily
            IAH1_BOVIN  (Q3SZ16)    , IAH1_DANRE  (Q503L4)    , IAH1_HUMAN  (Q2TAA2)    ,
            IAH1_MOUSE  (Q9DB29)    , IAH1_RAT    (Q711G3)    , IAH1_SCHPO  (O74648)    ,
            IAH1_YEAST  (P41734)    

        'GDSL' lipolytic enzyme family. Phospholipase B1 subfamily
            PLB1_CAVPO  (O70320)    , PLB1_HUMAN  (Q6P1J6)    , PLB1_MONDO  (Q06HQ7)    ,
            PLB1_MOUSE  (Q3TTY0)    , PLB1_RABIT  (Q05017)    , PLB1_RAT    (O54728)    

    Returns:
        [["'GDSL' lipolytic enzyme family", 'Q3MKY2'], ["'GDSL' lipolytic enzyme family", 'P86830'] ...]
    """
    uniprot_id_pat = re.compile(r'\([A-Z0-9]{4,}\)')
    raw_uniprot2protein_families = defaultdict(set)
    protein_family_name = ''
    with open(similar_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if 'family' in line:
                protein_family_name = line
            elif not line:
                protein_family_name = ''
            elif protein_family_name:
                uniprot_ids = uniprot_id_pat.findall(line)
                if uniprot_ids:
                    for _uniprot_id in uniprot_ids:
                        uniprot_id = _uniprot_id[1:-1]
                        raw_uniprot2protein_families[uniprot_id].add(protein_family_name)
    uniprot2protein_families = {}
    more_one_families_count = 0
    for k, v in raw_uniprot2protein_families.items():
        if len(v) > 1:
            more_one_families_count += 1
            logger.info(f'Uniprot_id {k} has >1 Protein_families {v}')
        uniprot2protein_families[k] = '; '.join(v)
    logger.info(f'len(merged_uniprot2protein_families) {len(uniprot2protein_families)}')
    uniprot2protein_families_df = pd.DataFrame(
        uniprot2protein_families.items(), columns =[UNIPROT_ID, PROTEIN_FAMILIES])
    if save_df:
        uniprot2protein_families_df.to_csv(uniprot2protein_families_file, index=False, sep=',')
        logger.info(f'more_one_families_count {more_one_families_count}')
        logger.info(uniprot2protein_families_df.shape)  # (511187, 2)
    return uniprot2protein_families, uniprot2protein_families_df


if __name__ == "__main__":

    # get_protein_families_with_uniport_ids('data/uniport/similar')
    # anaylyze_pdb_chain_uniprot()
    # parse_pdb_chain_uniprot()
    # parse_pdb_chain_uniprot_seq_fragments()
    # read_uniprot_sprot_fasta()
    # create_uniprot_to_seq_with_prot_families()
    # get_uniprot_seq_from_pdb_id_and_chain()
    get_pdb_id_and_chain_to_uniprot_id_dict()
    # parse_pdb_chain_uniprot()