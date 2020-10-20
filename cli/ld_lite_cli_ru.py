__version__ = 'V1.1'

from argparse import ArgumentParser, RawTextHelpFormatter

def add_args_ru(ver):
        '''
        Работа с аргументами командной строки.
        '''
        argparser = ArgumentParser(description=f'''
Программа выводит в табличном виде LD и дистанцию
между двумя вариантами, а также основополагающие
характеристики каждого из этих вариантов.

Версия: {ver}
Требуемые сторонние компоненты: pysam, tabulate
Автор: Платон Быкадоров (platon.work@gmail.com), 2020
Лицензия: GNU General Public License version 3
Поддержать проект: https://www.tinkoff.ru/rm/bykadorov.platon1/7tX2Y99140/
Документация: https://github.com/PlatonB/ld-tools/blob/master/README.md
Багрепорты, пожелания, общение: https://github.com/PlatonB/ld-tools/issues

Инструментарий ld_tools использует для
вычисления LD данные проекта 1000 Genomes.
Их скачивание и процессинг может занять много
времени, но осуществляется лишь однократно.

Результаты можно перенаправлять в
файл стандартным способом (> out.txt),
но целостность форматирования
сохраняемых таблиц не гарантирую.

Условные обозначения в справке по CLI:
- краткая форма с большой буквы: обязательный аргумент;
- в квадратных скобках: значение по умолчанию;
- в фигурных скобках: перечисление возможных значений.
''',
                                   formatter_class=RawTextHelpFormatter)
        argparser.add_argument('rs_id_1', metavar='str', type=str,
                               help='rsID 1-го варианта')
        argparser.add_argument('rs_id_2', metavar='str', type=str,
                               help='rsID 2-го варианта')
        argparser.add_argument('-D', '--intgen-dir-path', metavar='str', dest='intgen_dir_path', type=str,
                               help='Путь к папке для данных 1000G')
        argparser.add_argument('-f', '--skip-intgen-data-ver', dest='skip_intgen_data_ver', action='store_true',
                               help='Не проверять укомплектованность данных 1000G (сразу приступать к основным вычислениям)')
        argparser.add_argument('-g', '--gend-names', metavar='[both]', choices=['male', 'female', 'both'], default='both', dest='gend_names', type=str,
                               help='{male, female, both} Гендерная принадлежность сэмплов 1000G (для отбора определяющих LD генотипов)')
        argparser.add_argument('-e', '--pop-names', metavar='[all]', default='all', dest='pop_names', type=str,
                               help='''Популяционная принадлежность сэмплов 1000G (через запятую без пробела; для отбора определяющих LD генотипов;
https://www.internationalgenome.org/faq/which-populations-are-part-your-study/)''')
        args = argparser.parse_args()
        return args
