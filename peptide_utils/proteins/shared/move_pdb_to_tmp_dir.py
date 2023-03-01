import re, random
from pathlib import Path
import json
import pandas as pd
import numpy as np
from pandas import DataFrame
import os, sys, shutil


pdb_id = '1en2'
for file in Path('data/pdb/pdb_files').glob('*.pdb'):
    if file.stem == pdb_id:
        shutil.copyfile(file, (Path('data/pdb/pdb_files_tmp') / file.name))
        break
