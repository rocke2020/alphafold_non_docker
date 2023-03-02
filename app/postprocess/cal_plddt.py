import re, random
from pathlib import Path
import json
import pandas as pd
import numpy as np
from pandas import DataFrame
import os, sys
from icecream import ic
sys.path.append(os.path.abspath('.'))
ic.configureOutput(includeContext=True, argToStringFunction=lambda _: str(_))
from utils.log_util import logger


af_out_root = Path('/mnt/sdc/af_out/')
af_input_root = Path('/mnt/sdc/af_input')
ORDER = 'order'  # ranking order in 'ranking_debug.json'

# plddt_quality
HIGH = 'high'
WELL = 'well'
CAUTIOUS = 'cautious'
DISORDER = 'disorder'


def cal_plddt_quality(plddt_score):
    """ https://alphafold.ebi.ac.uk/faq
    
    Regions with pLDDT > 90 are expected to be modelled to high accuracy. These should be suitable for any application that benefits from high accuracy (e.g. characterising binding sites).
    Regions with pLDDT between 70 and 90 are expected to be modelled well (a generally good backbone prediction).
    Regions with pLDDT between 50 and 70 are low confidence and should be treated with caution.
    The 3D coordinates of regions with pLDDT < 50 often have a ribbon-like appearance and should not be interpreted. We show in our paper that pLDDT < 50 is a reasonably strong predictor of disorder, i.e. it suggests such a region is either unstructured in physiological conditions or only structured as part of a complex. (Note: this relationship has typically been tested in the context of well-studied proteins, which may have more evolutionarily-related sequences available than a randomly chosen UniProt entry.)
    Structured domains with many inter-residue contacts are likely to be more reliable than extended linkers or isolated long helices.
    Unphysical bond lengths and clashes do not usually appear in confident regions. Any part of a structure with several of these should be disregarded.
    Note that the PDB and mmCIF files contain coordinates for all regions, regardless of their pLDDT score. It is up to the user to interpret the model judiciously, in accordance with the guidance above.    
    """
    if plddt_score > 90:
        return HIGH
    elif plddt_score > 70:
        return WELL
    elif plddt_score >= 50:
        return CAUTIOUS
    else:
        return DISORDER


def cal_plddt(sub_dir_name='6itm_A_B_pdb_fasta_pos_human'):
    """

    """
    out_dir = af_out_root / sub_dir_name
    ranking_file = out_dir / 'ranking_debug.json'
    with open(ranking_file, 'r', encoding='utf-8') as f:
        ranking = json.load(f)

    inputs = read_multimer_inputs(filename_stem=sub_dir_name)
    prot_seq = inputs[0][1]
    pep_seq = inputs[1][1]
    ic(len(prot_seq))
    ic(len(pep_seq))
    best_model_name = ranking[ORDER][0]
    get_result_by_model_name(best_model_name, out_dir)


def read_multimer_inputs(filename_stem):
    """ multimer_input has 2 sequences in one fasta file
    prot first

    dummy example
    >6itm_A
    MGHHHHHHGSTELTPDQQTLLHFIMDSYNKQRMEITNKILKEAFSAEENFLILTEMATNHVQVLVEFTKKLPGFQTLDHEDQIALLKGSAVEAMFLRSVNDHKFTPLLCEIWDVQ
    >6itm_B
    KDHQLLRYLLDKDE

    """
    fasta_file = af_input_root / (filename_stem+'.fasta')
    # pdb_id_chain: pdb_fasta_seq
    inputs = []
    with open(fasta_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                pdb_id_chain = line[1:]
            elif line:
                inputs.append((pdb_id_chain, line))
    return inputs


def get_result_by_model_name(model_name, out_dir: Path):
    """
    the result_model_name.pkl is a dictionary, result_data.keys(): ['distogram', 'experimentally_resolved', 'masked_msa', 'num_recycles', 'predicted_aligned_error', 'predicted_lddt', 'structure_module', 'plddt', 'aligned_confidence_probs', 'max_predicted_aligned_error', 'ptm', 'iptm', 'ranking_confidence']

    num_recycles: 2, default

    predicted_lddt, dict, with only 1 key
        {'logits': DeviceArray([[-4.72423458e+00, -2.67970061e+00, -1.61319041e+00, ...,
                        -5.20858479e+00, -6.67728090e+00, -4.84406948e+00],
                        [-5.06980133e+00, -2.97154427e+00, -1.77395773e+00, ...,
                        -4.49883175e+00, -5.28827143e+00, -4.68708420e+00],
                        [-6.94888592e+00, -4.21163750e+00, -2.38282299e+00, ...,
                        -4.80423546e+00, -5.74000168e+00, -5.22399092e+00],
                        ...,
                        [-1.37081785e+01, -1.14596548e+01, -9.50779629e+00, ...,
                        1.84147787e+00,  1.99832344e+00,  2.14884996e+00],
                        [-1.10256052e+01, -1.09194975e+01, -9.36017799e+00, ...,
                        6.69315457e-03, -1.39155775e-01, -2.04017252e-01],
                        [-8.25554562e+00, -8.16820049e+00, -7.40738773e+00, ...,
                        -2.28720784e+00, -2.66509986e+00, -3.00252104e+00]],            dtype=float32)}

    plddt, np ndarray, 6itm_A_B_pdb_fasta
        [29.00279175 30.28257371 30.84997881 31.56364956 33.5245049  37.54542365
            39.26429317 42.33410803 50.44754624 62.20440454 78.26187201 89.01220113
            95.68556464 96.35764653 96.51805512 96.91498494 97.49222681 97.37942331
            97.13243092 97.61563768 97.87403913 97.11879553 96.22146629 97.44352593
            97.06477932 95.96894447 94.43797673 95.53276254 92.58950743 87.16342493
            77.14358865 75.09541624 75.58338486 79.49805735 81.33811475 84.40659475
            82.97362968 82.43925498 86.85816455 87.75583556 86.02434173 85.67132246
            88.04170497 86.90167663 88.36906858 91.23779824 93.07167343 90.90491956
            91.07698753 92.7225356  93.48038754 92.11299865 92.32953934 92.15849607
            91.72292475 92.88895511 93.63405128 91.73606683 93.78047114 95.53114146
            93.50726487 93.2696496  96.15096579 95.09669799 93.95769347 96.49970095
            96.74171788 95.24266536 96.45899042 97.44145262 95.48613567 94.28635761
            97.22535113 96.92666594 97.63489887 98.16293462 96.7234173  97.51182106
            98.19751781 97.34634881 95.48070323 96.37329581 98.09152648 98.31635792
            97.91787492 97.78502779 98.6082505  98.34059574 97.11627883 97.0438414
            97.64114115 97.47016063 96.55404738 97.45434192 97.16829954 95.06433938
            95.05118117 95.43836438 91.96308828 91.90155081 93.95411465 92.74580602
            87.88764973 89.8904387  90.7771878  83.09352689 73.97145375 71.0630423
            72.74271109 71.37792043 74.7943113  82.0084043  86.17783418 88.36550401
            88.91842294 90.15129708 91.667291   92.49288535 92.33281571 92.07755546
            92.41672859 92.3831825  90.13472073 86.28212618 86.7191053  86.09754685
            89.32525696 90.215315   89.38804371 89.98150285 92.7567071  92.4952756
            92.43039858 93.48942822 94.37596826 95.0802072  94.66086159 93.41385198
            94.33381306 94.59504676 93.78544244 93.73971312 94.08210779 92.19152865
            94.31674538 93.71418191 94.48946814 96.88984538 96.89492388 96.47904675
            97.19768672 98.30339679 97.98424991 97.90337908 98.62123176 98.59551022
            98.2330005  98.4943702  98.70228867 98.31683844 98.26545812 97.29434626
            96.70386376 93.35124194 95.80398064 97.91220942 97.91422938 98.43749993
            98.41115411 98.39579274 98.65327851 98.70325272 98.55249269 98.3417008
            98.4691596  98.56252384 98.1860171  98.1497513  98.41614931 97.93259873
            97.80533317 97.84948667 97.79122457 97.00952436 96.19426255 96.37158045
            95.41390618 92.6545059  91.69855705 89.31466034 89.46178885 85.41702446
            90.29532355 90.65480208 89.09769695 94.14344    96.00747139 96.49281126
            96.09922895 96.12195251 96.43017092 95.79275178 95.0424989  93.77238764
            94.39465707 92.26375301 92.87703681 93.8771821  93.71423208 90.38349161
            92.71564141 91.82364798 88.96802574 88.3327067  90.61172886 89.15772133
            87.28425226 85.4245224  87.29414288 84.68256557 82.17407641 83.7450894
            82.76731692 81.05127024 82.65594775 82.56180699 85.1306429  88.03575269
            92.32777917 95.28944458 95.60940114 97.18356918 96.20816567 94.53812701
            95.98918552 96.53134851 94.71597579 92.24348691 86.08075885 75.50365063
            60.62843742 78.59002393 92.8629017  95.03803393 96.67886669 96.99004461
            96.24528894 96.56067738 96.01454918 96.00950986 94.52092734 90.23315239
            78.72543128 61.9999746 ]
    """
    for file in out_dir.glob('*.pkl'):
        if file.stem.endswith(model_name):
            result_data = pd.read_pickle(file)
            ic(list(result_data.keys()))
            for k, v in result_data.items():
                if k == 'predicted_lddt':
                    predicted_lddt = v
                elif k == 'plddt':
                    plddt = v


if __name__ == "__main__":
    cal_plddt()
    pass