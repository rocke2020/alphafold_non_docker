import re, random
from pathlib import Path
import json
import os, sys, shutil


orig_out_dir = Path('/mnt/sdc/af_out')
simple_out_dir = Path('/mnt/sdc/af_out_simple')
simple_out_dir.mkdir(exist_ok=1)
existed_simple_sub_dirs = list(simple_out_dir.iterdir())


for file_dir in orig_out_dir.iterdir():
    new_sub_dir = simple_out_dir / file_dir.name
    new_sub_dir.mkdir(exist_ok=1)
    for file in file_dir.iterdir():
        if (file.stem.startswith('ranked_') and file.suffix == '.pdb') or file.stem == 'ranking_debug' or (
            file.stem == 'timings'
        ):
            new_file = new_sub_dir / file.name
            shutil.copyfile(file, new_file)

print('end')