__version__ = 'v1.0'

import os, shutil
from argparse import ArgumentParser, RawTextHelpFormatter

argparser = ArgumentParser(description=f'''
The script copies the results of
ld-tools (or other programs) scattered
in different subfolders into one folder.

Version: {__version__}
Dependencies: -
Author: Platon Bykadorov (platon.work@gmail.com), 2021
License: GNU General Public License version 3
Donate: https://www.tinkoff.ru/rm/bykadorov.platon1/7tX2Y99140/
Documentation: -
Bug reports, suggestions, talks: https://github.com/PlatonB/ld-tools/issues
''',
                           formatter_class=RawTextHelpFormatter)
argparser.add_argument('src_dir_path', metavar='str', type=str,
                       help='Path to folder with nested files')
argparser.add_argument('trg_dir_path', metavar='str', type=str,
                       help='Path to target folder')
args = argparser.parse_args()
for parent_dir_path, offspring_dir_names, offspring_file_names in os.walk(args.src_dir_path):
        if offspring_file_names == []:
                continue
        for offspring_file_name in offspring_file_names:
                shutil.copy2(os.path.join(parent_dir_path,
                                          offspring_file_name),
                             args.trg_dir_path)
