import re, random
from pathlib import Path
import json
import os, sys
sys.path.append(os.path.abspath('.'))
from icecream import ic
ic.configureOutput(includeContext=True, argToStringFunction=lambda _: str(_))
from utils.file_util import FileUtil


def add_2_blanks_at_the_end_for_md_file(file_path):
    """  """
    lines = FileUtil.read_raw_text(file_path)
    lines_with_2_blanks_at_the_end = []
    for line in lines:
        if not line.startswith('#') and line != '':
            line = line + '  '
        lines_with_2_blanks_at_the_end.append(line)
    FileUtil.write_raw_text(lines_with_2_blanks_at_the_end, file_path)


if __name__ == "__main__":
    add_2_blanks_at_the_end_for_md_file('data_process/tasks/anti_inflammatory/readme.md')
    pass