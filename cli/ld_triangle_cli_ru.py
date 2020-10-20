__version__ = 'V1.1'

from argparse import ArgumentParser, RawTextHelpFormatter

def add_args_ru(ver):
        '''
        Работа с аргументами командной строки.
        '''
        argparser = ArgumentParser(description=f'''
Программа, строящая LD-матрицы для
всех пар каждого набора вариантов в виде
треугольной тепловой карты и/или таблицы.

Версия: {ver}
Требуемые сторонние компоненты: pysam, plotly
Автор: Платон Быкадоров (platon.work@gmail.com), 2018-2020
Лицензия: GNU General Public License version 3
Поддержать проект: https://www.tinkoff.ru/rm/bykadorov.platon1/7tX2Y99140/
Документация: https://github.com/PlatonB/ld-tools/blob/master/README.md
Багрепорты, пожелания, общение: https://github.com/PlatonB/ld-tools/issues

Поддерживаемые исходные файлы - таблицы,
содержащие столбец с набором rsIDs. Если
таких столбцов - несколько, программа
будет использовать самый левый.

В одном исходном файле могут
содержаться данные разных хромосом.
Для каждой хромосомы программа
построит отдельную матрицу.

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
                               help='Путь к папке для данных 1000G')
        argparser.add_argument('-t', '--trg-top-dir-path', metavar='[None]', dest='trg_top_dir_path', type=str,
                               help='Путь к папке для результатов (по умолчанию - путь к исходной папке)')
        argparser.add_argument('-m', '--meta-lines-quan', metavar='[0]', default=0, dest='meta_lines_quan', type=int,
                               help='Количество строк метаинформации, включая шапку, в начале каждой исходной таблицы')
        argparser.add_argument('-f', '--skip-intgen-data-ver', dest='skip_intgen_data_ver', action='store_true',
                               help='Не проверять укомплектованность данных 1000G (сразу приступать к основным вычислениям)')
        argparser.add_argument('-g', '--gend-names', metavar='[both]', choices=['male', 'female', 'both'], default='both', dest='gend_names', type=str,
                               help='{male, female, both} Гендерная принадлежность сэмплов 1000G (для отбора определяющих LD генотипов)')
        argparser.add_argument('-e', '--pop-names', metavar='[all]', default='all', dest='pop_names', type=str,
                               help='''Популяционная принадлежность сэмплов 1000G (через запятую без пробела; для отбора определяющих LD генотипов;
https://www.internationalgenome.org/faq/which-populations-are-part-your-study/)''')
        argparser.add_argument('-l', '--ld-measure', metavar='[r_square]', choices=['r_square', 'd_prime'], default='r_square', dest='ld_measure', type=str,
                               help='{r_square, d_prime} Мера LD для построения матриц и выставления нижнего порога LD')
        argparser.add_argument('-z', '--ld-low-thres', metavar='[None]', dest='ld_low_thres', type=float,
                               help='Нижний порог LD (подпороговые значения приравняются к нулю)')
        argparser.add_argument('-o', '--matrix-type', metavar='[heatmap]', choices=['heatmap', 'table', 'both'], default='heatmap', dest='matrix_type', type=str,
                               help='{heatmap, table, both} Вид матриц значений LD')
        argparser.add_argument('-j', '--heatmap-json', dest='heatmap_json', action='store_true',
                               help='Сохранение объектов тепловых карт в виде JSON (полезно для дебага)')
        argparser.add_argument('-i', '--disp-letters', dest='disp_letters', action='store_true',
                               help='Вывод LD-значений и rsID-лейблов осей на тепловую карту')
        argparser.add_argument('-c', '--color-pal', metavar='[greens]', default='greens', dest='color_pal', type=str,
                               help='Цветовая палитра тепловой карты (https://github.com/PlatonB/ld-tools#подходящие-цветовые-палитры)')
        argparser.add_argument('-k', '--font-size', metavar='[None]', dest='font_size', type=int,
                               help='Размер шрифта надписей на тепловой карте (plotly: по умолчанию - 12; больше диаграмма - мельче делайте шрифт)')
        argparser.add_argument('-q', '--square-shape', dest='square_shape', action='store_true',
                               help='Квадратная форма тепловой карты')
        argparser.add_argument('-s', '--dont-disp-footer', dest='dont_disp_footer', action='store_true',
                               help='Не выводить на тепловую карту информацию о программе')
        argparser.add_argument('-p', '--max-proc-quan', metavar='[4]', default=4, dest='max_proc_quan', type=int,
                               help='Максимальное количество параллельно обрабатываемых таблиц')
        args = argparser.parse_args()
        return args
