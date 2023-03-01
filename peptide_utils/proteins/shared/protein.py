import dataclasses
import io
from typing import Any, Mapping, Optional
import peptide_utils.proteins.shared.residue_constants as residue_constants
from Bio.PDB import PDBParser
import numpy as np
from string import ascii_uppercase, ascii_lowercase
alphabet_list = list(ascii_uppercase+ascii_lowercase)
from icecream import ic
ic.configureOutput(includeContext=True, argToStringFunction=lambda _: str(_))
from utils.log_util import logger


MODRES = {'MSE':'MET','MLY':'LYS','FME':'MET','HYP':'PRO',
          'TPO':'THR','CSO':'CYS','SEP':'SER','M3L':'LYS',
          'HSK':'HIS','SAC':'SER','PCA':'GLU','DAL':'ALA',
          'CME':'CYS','CSD':'CYS','OCS':'CYS','DPR':'PRO',
          'B3K':'LYS','ALY':'LYS','YCM':'CYS','MLZ':'LYS',
          '4BF':'TYR','KCX':'LYS','B3E':'GLU','B3D':'ASP',
          'HZP':'PRO','CSX':'CYS','BAL':'ALA','HIC':'HIS',
          'DBZ':'ALA','DCY':'CYS','DVA':'VAL','NLE':'LEU',
          'SMC':'CYS','AGM':'ARG','B3A':'ALA','DAS':'ASP',
          'DLY':'LYS','DSN':'SER','DTH':'THR','GL3':'GLY',
          'HY3':'PRO','LLP':'LYS','MGN':'GLN','MHS':'HIS',
          'TRQ':'TRP','B3Y':'TYR','PHI':'PHE','PTR':'TYR',
          'TYS':'TYR','IAS':'ASP','GPL':'LYS','KYN':'TRP',
          'CSD':'CYS','SEC':'CYS'}

def pdb_to_string(pdb_file):
  modres = {**MODRES}
  lines = []
  seen = []
  for line in open(pdb_file,"rb"):
    line = line.decode("utf-8","ignore").rstrip()
    if line[:6] == "MODRES":
      k = line[12:15]
      v = line[24:27]
      if k not in modres and v in residue_constants.restype_3to1:
        modres[k] = v
    if line[:6] == "HETATM":
      k = line[17:20]
      if k in modres:
        line = "ATOM  "+line[6:17]+modres[k]+line[20:]
    if line[:4] == "ATOM":
      chain = line[21:22]
      atom = line[12:12+4].strip()
      resi = line[17:17+3]
      resn = line[22:22+5].strip()
      if resn[-1].isalpha(): # alternative atom
        resn = resn[:-1]
        line = line[:26]+" "+line[27:]
      key = f"{chain}_{resn}_{resi}_{atom}"
      if key not in seen: # skip alternative placements
        lines.append(line)
        seen.append(key)
  return "\n".join(lines)


@dataclasses.dataclass(frozen=True)
class Protein:
  """Protein structure representation."""

  # Cartesian coordinates of atoms in angstroms. The atom types correspond to
  # residue_constants.atom_types, i.e. the first three are N, CA, CB.
  atom_positions: np.ndarray  # [num_res, num_atom_type, 3]

  # Amino-acid type for each residue represented as an integer between 0 and
  # 20, where 20 is 'X'.
  aatype: np.ndarray  # [num_res]

  # Binary float mask to indicate presence of a particular atom. 1.0 if an atom
  # is present and 0.0 if not. This should be used for loss masking.
  atom_mask: np.ndarray  # [num_res, num_atom_type]

  # Residue index as used in PDB. It is not necessarily continuous or 0-indexed.
  residue_index: np.ndarray  # [num_res]

  # B-factors, or temperature factors, of each residue (in sq. angstroms units),
  # representing the displacement of the residue from its ground truth mean
  # value.
  b_factors: np.ndarray  # [num_res, num_atom_type]


def from_pdb_file(pdb_fh: str, chain_id: Optional[str] = None) -> Protein:
  """Takes a PDB string and constructs a Protein object.

  WARNING: All non-standard residue types will be converted into UNK. All
    non-standard atoms will be ignored.

  Args:
    pdb_str: The contents of the pdb file
    chain_id: If None, then the pdb file must contain a single chain (which
      will be parsed). If chain_id is specified (e.g. A), then only that chain
      is parsed.

  Returns:
    A new `Protein` parsed from the pdb contents.
  """
  parser = PDBParser(QUIET=True)
  structure = parser.get_structure('none', pdb_fh)
  models = list(structure.get_models())
  model = models[0]
  if chain_id is not None:
    chain = model[chain_id]
  else:
    chains = list(model.get_chains())
    if len(chains) != 1:
      raise ValueError(
          'Only single chain PDBs are supported when chain_id not specified. '
          f'Found {len(chains)} chains.')
    else:
      chain = chains[0]

  atom_positions = []
  aatype = []
  atom_mask = []
  residue_index = []
  b_factors = []

  # only choose the last disorded res for this loop.
  for res in chain:
    # ic(res)
    hetflag, resseq, icode = res.get_id()
    if hetflag != ' ': continue
    res_shortname = residue_constants.restype_3to1.get(res.resname, 'X')
    restype_idx = residue_constants.restype_order.get(
        res_shortname, residue_constants.restype_num)
    pos = np.zeros((residue_constants.atom_type_num, 3))
    mask = np.zeros((residue_constants.atom_type_num,))
    res_b_factors = np.zeros((residue_constants.atom_type_num,))
    for atom in res:
      if atom.name not in residue_constants.atom_types:
        continue
      pos[residue_constants.atom_order[atom.name]] = atom.coord
      mask[residue_constants.atom_order[atom.name]] = 1.
      res_b_factors[residue_constants.atom_order[atom.name]] = atom.bfactor
    if np.sum(mask) < 0.5:
      # If no known atom positions are reported for the residue then skip it.
      continue
    aatype.append(restype_idx)
    atom_positions.append(pos)
    atom_mask.append(mask)
    residue_index.append(res.id[1])
    b_factors.append(res_b_factors)

  return Protein(
      atom_positions=np.array(atom_positions),
      atom_mask=np.array(atom_mask),
      aatype=np.array(aatype),
      residue_index=np.array(residue_index),
      b_factors=np.array(b_factors))


def get_seq_from_pdb_file_pdb_parser(pdb_fh: str, chain_id: Optional[str] = None) -> Protein:
  """Takes a PDB string and constructs a Protein object.

  WARNING: All non-standard residue types will be converted into UNK. All non-standard atoms will be ignored.

  Args:
    pdb_str: The contents of the pdb file
    chain_id: If None, then the pdb file must contain a single chain (which
      will be parsed). If chain_id is specified (e.g. A), then only that chain
      is parsed.

  Returns:
    A aa sequence which only indluce the disordered aa.
  """
  # logger.info(f'pdb_id {pdb_fh} {chain_id}')
  parser = PDBParser(QUIET=True)
  structure = parser.get_structure('none', pdb_fh)
  models = list(structure.get_models())
  model = models[0]
  if chain_id is not None:
    chain = model[chain_id]
  else:
    chains = list(model.get_chains())
    if len(chains) != 1:
      raise ValueError(
          'Only single chain PDBs are supported when chain_id not specified. '
          f'Found {len(chains)} chains.')
    else:
      chain = chains[0]

  aas = chain.get_unpacked_list()
  aas_no_disordered = []
  _unique_id = ''
  for res in aas:
    # <Residue ILE het=  resseq=16 icode= >, <Residue VAL het=  resseq=17 icode= > 
    ## disorder residue because duplicate f'{resseq}_{icode}'
    # <Residue CYS het=  resseq=22 icode= >, <Residue SER het=  resseq=22 icode= >
    ## het
    # <Residue CA het=H_CA resseq=1246 icode= >, <Residue HOH het=W resseq=2001 icode= >
    hetflag, resseq, icode = res.get_id()
    if hetflag != ' ':
      # ic(res)
      continue
    res_shortname = residue_constants.restype_3to1.get(res.resname, 'X')
    unique_id = f'{resseq}_{icode}'
    ## disorder residue because duplicate f'{resseq}_{icode}', only keep the first one to align with FASTA format
    if _unique_id == unique_id:
      continue
    _unique_id = unique_id
    # ic(res_shortname)
    aas_no_disordered.append(res_shortname)
  return ''.join(aas_no_disordered)


def renum_pdb_str(pdb_str, Ls=None, renum=True, offset=1):
  if Ls is not None:
    L_init = 0
    new_chain = {}
    for L,c in zip(Ls, alphabet_list):
      new_chain.update({i:c for i in range(L_init,L_init+L)})
      L_init += L  

  n,num,pdb_out = 0,offset,[]
  resnum_ = None
  chain_ = None
  new_chain_ = new_chain[0]
  for line in pdb_str.split("\n"):
    if line[:4] == "ATOM":
      chain = line[21:22]
      resnum = int(line[22:22+5])
      if resnum_ is None: resnum_ = resnum
      if chain_ is None: chain_ = chain
      if resnum != resnum_ or chain != chain_:
        num += (resnum - resnum_)  
        n += 1
        resnum_,chain_ = resnum,chain
      if Ls is not None:
        if new_chain[n] != new_chain_:
          num = offset
          new_chain_ = new_chain[n]
      N = num if renum else resnum
      if Ls is None: pdb_out.append("%s%4i%s" % (line[:22],N,line[26:]))
      else: pdb_out.append("%s%s%4i%s" % (line[:21],new_chain[n],N,line[26:]))        
  return "\n".join(pdb_out)
