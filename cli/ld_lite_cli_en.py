__version__ = 'V1.0'

from argparse import ArgumentParser, RawTextHelpFormatter

def add_args_en(ver):
        '''
        Работа с аргументами командной строки.
        '''
        argparser = ArgumentParser(description=f'''
The program prints in tabular form LD and the
distance between the two variants, as well as the
essential characteristics of each of these variants.

Version: {ver}
Dependencies: pysam, tabulate
Author: Platon Bykadorov (platon.work@gmail.com), 2020
License: GNU General Public License version 3
Donate: https://www.tinkoff.ru/rm/bykadorov.platon1/7tX2Y99140/
Documentation: https://github.com/PlatonB/ld-tools/blob/master/README-EN.md
Bug reports, suggestions, talks: https://github.com/PlatonB/ld-tools/issues

ld_tools uses 1000 Genomes project data for LD
calculation. Their downloading and processing
may take a long time, but it is done only once.

The results can be redirected to a
file in a standard way (> out.txt),
but the integrity of the formatting
of the saved tables is not guaranteed.

CLI help legend:
- a short form with a capital letter: mandatory argument;
- in square brackets: default value;
- in curly brackets: list of possible values.
''',
                                   formatter_class=RawTextHelpFormatter)
        argparser.add_argument('rs_id_1', metavar='str', type=str,
                               help='rsID of the first variant')
        argparser.add_argument('rs_id_2', metavar='str', type=str,
                               help='rsID of the second variant')
        argparser.add_argument('-D', '--intgen-dir-path', metavar='str', dest='intgen_dir_path', type=str,
                               help='Path to folder for 1000 Genomes data')
        argparser.add_argument('-f', '--skip-intgen-data-ver', dest='skip_intgen_data_ver', action='store_true',
                               help='Do not check 1000 Genomes data completeness (start main calculations immediately)')
        argparser.add_argument('-g', '--gend-names', metavar='[both]', choices=['male', 'female', 'both'], default='both', dest='gend_names', type=str,
                               help='{male, female, both} Belonging of 1000 Genomes samples to genders')
        argparser.add_argument('-e', '--pop-names', metavar='[all]', default='all', dest='pop_names', type=str,
                               help='Belonging of 1000 Genomes samples to populations (separated by commas without space; https://www.internationalgenome.org/faq/which-populations-are-part-your-study/)')
        args = argparser.parse_args()
        return args
