__version__ = 'V1.0'

from argparse import ArgumentParser, RawTextHelpFormatter

def add_args_ru(ver):
        '''
        Работа с аргументами командной строки.
        '''
        argparser = ArgumentParser(description=f'''
Программа ищет в пределах фланков варианты,
обладающие надпороговым значением неравновесия
по сцеплению с каждым запрашиваемым вариантом.

Версия: {ver}
Требуемый сторонний компонент: pysam
Автор: Платон Быкадоров (platon.work@gmail.com), 2018-2020
Лицензия: GNU General Public License version 3
Поддержать проект: https://www.tinkoff.ru/rm/bykadorov.platon1/7tX2Y99140/
Документация: https://github.com/PlatonB/ld-tools/blob/master/README.md
Багрепорты, пожелания, общение: https://github.com/PlatonB/ld-tools/issues

Поддерживаемые исходные файлы - таблицы,
содержащие столбец с набором rsIDs. Если
таких столбцов - несколько, программа
будет использовать самый левый.

Инструментарий ld_tools использует для
вычисления LD данные проекта 1000 Genomes.
Их скачивание и процессинг может занять много
времени, но осуществляется лишь однократно.

Условные обозначения в справке по CLI:
- краткая форма с большой буквы: обязательный аргумент;
- в квадратных скобках: значение по умолчанию;
- в фигурных скобках: перечисление возможных значений.
''',
                                   formatter_class=RawTextHelpFormatter)
        argparser.add_argument('-S', '--src-dir-path', metavar='str', dest='src_dir_path', type=str,
                               help='Путь к папке с исходными таблицами')
        argparser.add_argument('-D', '--intgen-dir-path', metavar='str', dest='intgen_dir_path', type=str,
                               help='Путь к папке для данных 1000 Genomes')
        argparser.add_argument('-t', '--trg-top-dir-path', metavar='[None]', dest='trg_top_dir_path', type=str,
                               help='Путь к папке для результатов (по умолчанию - путь к исходной папке)')
        argparser.add_argument('-m', '--meta-lines-quan', metavar='[0]', default=0, dest='meta_lines_quan', type=int,
                               help='Количество строк метаинформации, включая шапку, в начале каждой исходной таблицы')
        argparser.add_argument('-f', '--skip-intgen-data-ver', dest='skip_intgen_data_ver', action='store_true',
                               help='Не проверять укомплектованность данных 1000 Genomes (сразу приступать к основным вычислениям)')
        argparser.add_argument('-g', '--gend-names', metavar='[both]', choices=['male', 'female', 'both'], default='both', dest='gend_names', type=str,
                               help='{male, female, both} Гендерная принадлежность сэмплов 1000 Genomes')
        argparser.add_argument('-e', '--pop-names', metavar='[all]', default='all', dest='pop_names', type=str,
                               help='Популяционная принадлежность сэмплов 1000 Genomes (через запятую без пробела; https://www.internationalgenome.org/faq/which-populations-are-part-your-study/)')
        argparser.add_argument('-w', '--flank-size', metavar='[100000]', default=100000, dest='flank_size', type=int,
                               help='Размер *каждого* из фланков, в пределах которых надо искать неравновесно сцепленные варианты')
        argparser.add_argument('-l', '--ld-thres-measure', metavar='[r_square]', choices=['r_square', 'd_prime'], default='r_square', dest='ld_thres_measure', type=str,
                               help='{r_square, d_prime} Мера для выставления нижнего порога LD')
        argparser.add_argument('-z', '--ld-low-thres', metavar='[0.8]', default=0.8, dest='ld_low_thres', type=float,
                               help='Нижний порог LD')
        argparser.add_argument('-o', '--trg-file-type', metavar='[tsv]', choices=['tsv', 'json', 'rsids'], default='tsv', dest='trg_file_type', type=str,
                               help='{tsv, json, rsids} Формат конечных файлов')
        argparser.add_argument('-p', '--max-proc-quan', metavar='[4]', default=4, dest='max_proc_quan', type=int,
                               help='Максимальное количество параллельно обрабатываемых таблиц')
        args = argparser.parse_args()
        return args
