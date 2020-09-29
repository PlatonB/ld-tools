__version__ = 'V1.0'

from argparse import ArgumentParser, RawTextHelpFormatter

def add_args_en(ver):
        '''
        Работа с аргументами командной строки.
        '''
        argparser = ArgumentParser(description=f'''
The program searches variants for each source
variant within the window and with a linkage
disequilibrium above the threshold value.

Version: {ver}
Dependency: pysam
Author: Platon Bykadorov (platon.work@gmail.com), 2018-2020
License: GNU General Public License version 3
Donate: https://www.tinkoff.ru/rm/bykadorov.platon1/7tX2Y99140/
Documentation: https://github.com/PlatonB/ld-tools/blob/master/README-EN.md
Bug reports, suggestions, talks: https://github.com/PlatonB/ld-tools/issues

Supported source files are tables containing
a column with rsIDs. If there are more than 1
rsID columns, the program will use the left one.

ld_tools uses 1000 Genomes project data for LD
calculation. Their downloading and processing
may take a long time, but it is done only once.

CLI help legend:
- a short form with a capital letter: mandatory argument;
- in square brackets: default value;
- in curly brackets: list of possible values.
''',
                                   formatter_class=RawTextHelpFormatter)
        argparser.add_argument('-S', '--src-dir-path', metavar='str', dest='src_dir_path', type=str,
                               help='Path to folder with source tables')
        argparser.add_argument('-D', '--intgen-dir-path', metavar='str', dest='intgen_dir_path', type=str,
                               help='Path to folder for 1000 Genomes data')
        argparser.add_argument('-t', '--trg-top-dir-path', metavar='[None]', dest='trg_top_dir_path', type=str,
                               help='Path to target folder (default: path to source folder)')
        argparser.add_argument('-m', '--meta-lines-quan', metavar='[0]', default=0, dest='meta_lines_quan', type=int,
                               help='Number of meta-information lines (including line with column names)')
        argparser.add_argument('-f', '--skip-intgen-data-ver', dest='skip_intgen_data_ver', action='store_true',
                               help='Do not check 1000 Genomes data completeness (start main calculations immediately)')
        argparser.add_argument('-g', '--gend-names', metavar='[both]', choices=['male', 'female', 'both'], default='both', dest='gend_names', type=str,
                               help='{male, female, both} Belonging of 1000 Genomes samples to genders')
        argparser.add_argument('-e', '--pop-names', metavar='[all]', default='all', dest='pop_names', type=str,
                               help='Belonging of 1000 Genomes samples to populations (separated by commas without space; https://www.internationalgenome.org/faq/which-populations-are-part-your-study/)')
        argparser.add_argument('-w', '--flank-size', metavar='[100000]', default=100000, dest='flank_size', type=int,
                               help='The size of each of the flanks, where to look for in-LD variants')
        argparser.add_argument('-l', '--ld-thres-measure', metavar='[r_square]', choices=['r_square', 'd_prime'], default='r_square', dest='ld_thres_measure', type=str,
                               help='{r_square, d_prime} Measure for setting the lower LD threshold')
        argparser.add_argument('-z', '--ld-low-thres', metavar='[0.8]', default=0.8, dest='ld_low_thres', type=float,
                               help='Lower LD threshold')
        argparser.add_argument('-o', '--trg-file-type', metavar='[tsv]', choices=['tsv', 'json', 'rsids'], default='tsv', dest='trg_file_type', type=str,
                               help='{tsv, json, rsids} Target file format')
        argparser.add_argument('-p', '--max-proc-quan', metavar='[4]', default=4, dest='max_proc_quan', type=int,
                               help='Maximum number of tables to be processed in parallel')
        args = argparser.parse_args()
        return args
