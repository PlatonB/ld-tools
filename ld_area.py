__version__ = 'V10.1'

def add_args():
        '''
        Работа с аргументами командной строки.
        '''
        argparser = ArgumentParser(description=f'''
Программа ищет в пределах фланков SNPs,
обладающие надпороговым значением неравновесия
по сцеплению с каждым запрашиваемым SNP.

Автор: Платон Быкадоров (platon.work@gmail.com), 2018-2020
Версия: {__version__}
Лицензия: GNU General Public License version 3
Поддержать проект: https://money.yandex.ru/to/41001832285976
Документация: https://github.com/PlatonB/ld-tools/blob/master/README.md
Багрепорты/пожелания/общение: https://github.com/PlatonB/ld-tools/issues

Перед запуском программы нужно
установить pysam (см. документацию).

Поддерживаемые исходные файлы - таблицы,
содержащие столбец с набором refSNPIDs.
Если таких столбцов - несколько,
программа будет использовать самый левый.

Инструментарий ld_tools использет для
вычисления LD данные проекта 1000 Genomes.
Их скачивание и процессинг может занять много
времени, но осуществляется лишь однократно.

Условные обозначения в справке по CLI:
- краткая форма с большой буквы - обязательный аргумент;
- в квадратных скобках - значение по умолчанию;
- в фигурных скобках - перечисление возможных значений.
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
                               help='Размер *каждого* из фланков, в пределах которых надо искать сцепленные SNP')
        argparser.add_argument('-l', '--ld-thres-measure', metavar='[r_square]', choices=['r_square', 'd_prime'], default='r_square', dest='ld_thres_measure', type=str,
                               help='{r_square, d_prime} Мера для выставления нижнего порога LD')
        argparser.add_argument('-z', '--ld-low-thres', metavar='[0.5]', default=0.5, dest='ld_low_thres', type=float,
                               help='Нижний порог LD')
        argparser.add_argument('-o', '--trg-file-type', metavar='[tsv]', choices=['tsv', 'json', 'rsids'], default='tsv', dest='trg_file_type', type=str,
                               help='{tsv, json, rsids} Формат конечных файлов')
        argparser.add_argument('-p', '--max-proc-quan', metavar='[4]', default=4, dest='max_proc_quan', type=int,
                               help='Максимальное количество параллельно обрабатываемых таблиц')
        args = argparser.parse_args()
        return args

def build_ucsc_header(header_key, header_val):
        '''
        Конечные не-JSON-файлы будут начинаться
        с хэдеров, структура которых напоминает
        таковую у хэдеров UCSC Table Browser.
        Эта функция собирает каждый элемент хэдера.
        '''
        if type(header_val).__name__ == 'str':
                header_val = f'"{header_val}"'
        elif type(header_val).__name__ == 'tuple':
                header_val = ','.join([f'"{header_val_element}"' for header_val_element in header_val])
        return f'{header_key}={header_val}'

class PrepSingleProc():
        '''
        Класс, спроектированный
        под безопасный параллельный
        поиск SNPs, сцепленных со
        снипами из исходного файла.
        '''
        def __init__(self, args):
                '''
                Получение атрибутов, необходимых
                заточенной под многопроцессовое
                выполнение функции получения SNPs,
                обладающих надпороговым LD по
                отношению к SNPs исходного набора.
                Атрибуты должны быть созданы
                единожды и далее ни в
                коем случае не изменяться.
                Получаются они в основном из
                указанных исследователем опций.
                '''
                self.src_dir_path = os.path.normpath(args.src_dir_path)
                self.intgen_dir_path = os.path.normpath(args.intgen_dir_path)
                if args.trg_top_dir_path is None:
                        self.trg_top_dir_path = self.src_dir_path
                else:
                        self.trg_top_dir_path = os.path.normpath(args.trg_top_dir_path)
                self.meta_lines_quan = args.meta_lines_quan
                if args.skip_intgen_data_ver:
                        self.intgen_convdb_path = os.path.join(self.intgen_dir_path, 'conversion.db')
                else:
                        self.intgen_convdb_path = prep_intgen_data(self.intgen_dir_path)
                if args.gend_names == 'male':
                        self.gend_names = ('male',)
                elif args.gend_names == 'female':
                        self.gend_names = ('female',)
                else:
                        self.gend_names = ('male', 'female')
                self.pop_names = tuple(args.pop_names.upper().split(','))
                self.sample_names = get_sample_names(self.gend_names,
                                                     self.pop_names,
                                                     self.intgen_convdb_path)
                self.flank_size = args.flank_size
                self.ld_thres_measure = args.ld_thres_measure
                self.ld_low_thres = args.ld_low_thres
                self.trg_file_type = args.trg_file_type
                
        def get_lnkd_snps(self, src_file_name):
                '''
                Функция нахождения SNPs,
                находящихся в неравновесии
                по сцеплению с каждым SNP.
                '''
                
                #Подготовка к считыванию основной
                #части файла исследователя.
                with open(os.path.join(self.src_dir_path, src_file_name)) as src_file_opened:
                        for meta_line_index in range(self.meta_lines_quan):
                                src_file_opened.readline()
                                
                        #Значительная часть программы
                        #привязана к функциональности pysam.
                        #Эта библиотека поддерживает запросы
                        #по геномным координатам, поэтому
                        #далее будем оперировать только ими.
                        #Извлечём имена хромосом и позиции
                        #исходных SNPs из конвертационной
                        #БД и накопим их в словарь, сделав
                        #последнее таким образом, чтобы первые
                        #стали ключами, а вторые - значениями.
                        #По ходу работы с базой данных программа
                        #скипнет невалидные исходные снипы.
                        pos_by_chrs = {}
                        with sqlite3.connect(self.intgen_convdb_path) as conn:
                                cursor = conn.cursor()
                                for line in src_file_opened:
                                        try:
                                                rs_id = re.search(r'rs\d+', line).group()
                                        except AttributeError:
                                                continue
                                        cursor.execute(f'SELECT CHROM, POS FROM snps WHERE ID IN ("{rs_id}")')
                                        chrom_pos = cursor.fetchone()
                                        if chrom_pos is None:
                                                continue
                                        try:
                                                pos_by_chrs[chrom_pos[0]].append(chrom_pos[1])
                                        except KeyError:
                                                pos_by_chrs[chrom_pos[0]] = [chrom_pos[1]]
                                cursor.close()
                                
                #Дальнейшая обработка данных текущего
                #файла возможна, только если в их составе
                #обнаружился хоть один валидный rsID.
                if len(pos_by_chrs) > 0:
                        
                        #В одну папку второго уровня
                        #планируется размещать все
                        #результаты, полученные по
                        #данным одного исходного файла.
                        src_file_base = '.'.join(src_file_name.split('.')[:-1])
                        trg_dir_path = os.path.join(self.trg_top_dir_path, f'{src_file_base}_lnkd')
                        os.mkdir(trg_dir_path)
                        
                        #Сокращённые заглавия таких
                        #глобальных характеристик,
                        #как хромосома и настройки
                        #отбора оппонирующих SNPs.
                        meta_keys = ['chr',
                                     'gends',
                                     'pops',
                                     'each_flank',
                                     f'{self.ld_thres_measure}_thres']
                        
                        #Имена столбцов TSV.
                        #Они же - ключи JSON.
                        header_row = ['hg38_pos',
                                      'rsID',
                                      'ref',
                                      'alt',
                                      'type',
                                      'alt_freq',
                                      'r2',
                                      "D'",
                                      'dist']
                        
                        #Далее работа будет проводиться
                        #для данных каждой хромосомы
                        #по-отдельности, а результаты
                        #пойдут в соответствующие хромосомные
                        #подпапки (третий уровень вложенности).
                        for chrom in pos_by_chrs:
                                chr_dir_path = os.path.join(trg_dir_path, chrom)
                                os.mkdir(chr_dir_path)
                                if self.trg_file_type in ['tsv', 'json']:
                                        ext = self.trg_file_type
                                else:
                                        ext = 'txt'
                                        
                                #Собственно, глобальные
                                #характеристики: текущая
                                #хромосома и заданные
                                #исследователем опции.
                                #Пойдут в TSV или JSON.
                                meta_vals = [chrom,
                                             self.gend_names,
                                             self.pop_names,
                                             self.flank_size,
                                             self.ld_low_thres]
                                
                                #Хэдер для TSV или файла с чистыми rsIDs.
                                #Его структура напоминает таковую у
                                #таблиц, генерируемых UCSC Table Browser.
                                ucsc_header_line = '##' + ' '.join(map(build_ucsc_header,
                                                                       meta_keys,
                                                                       meta_vals))
                                
                                #Поскольку известно, какая хромосома
                                #сейчас обрабатывается, можно сразу
                                #открыть именно относящимися к этой
                                #хромосоме архив с 1000 Genomes-данными.
                                #Из него для начала будем извлекать
                                #аннотации только запрашиваемых SNPs.
                                with VariantFile(os.path.join(self.intgen_dir_path, f'{chrom}.vcf.gz')) as intgen_vcf_opened:
                                        for pos in pos_by_chrs[chrom]:
                                                for query_snp_data in intgen_vcf_opened.fetch(chrom, pos - 1, pos):
                                                        trg_file_name = f'{chrom}_{query_snp_data.id}_{self.ld_thres_measure[0]}_{str(self.ld_low_thres)}.{ext}'
                                                        trg_file_path = os.path.join(chr_dir_path, trg_file_name)
                                                        if os.path.exists(trg_file_path):
                                                                continue
                                                        
                                                        #Создаём флаг, по которому далее будет
                                                        #определено, оказались ли в конечном
                                                        #файле строки, отличные от хэдеров.
                                                        empty_res = True
                                                        
                                                        #Вычисляем границы фланков вокруг
                                                        #текущего запрашиваемого SNP, где
                                                        #надо искать сцепленные с ним SNPs.
                                                        #Если нижняя граница заползёт за
                                                        #ноль, то во избежание ошибки со
                                                        #стороны pysam приравняем её к нулю.
                                                        low_bound = query_snp_data.pos - self.flank_size
                                                        if low_bound < 0:
                                                                low_bound = 0
                                                        high_bound = query_snp_data.pos + self.flank_size
                                                        
                                                        #Получение генотипов и
                                                        #ключевых характеристик
                                                        #запрашиваемого SNP.
                                                        query_snp_genotypes = []
                                                        for sample_name in self.sample_names:
                                                                try:
                                                                        query_snp_genotypes += query_snp_data.samples[sample_name]['GT']
                                                                except KeyError:
                                                                        continue
                                                        query_snp_alt_freq = round(query_snp_genotypes.count(1) /
                                                                                   len(query_snp_genotypes), 4)
                                                        query_snp_ann = [query_snp_data.pos,
                                                                         query_snp_data.id,
                                                                         query_snp_data.ref,
                                                                         ','.join(query_snp_data.alts),
                                                                         ','.join(query_snp_data.info['VT']),
                                                                         query_snp_alt_freq] + ['quer'] * 3
                                                        
                                                        #Открываем конечный файл на запись
                                                        #и проделываем кропотливую работу по
                                                        #созданию формат-ориентированных хэдеров.
                                                        with open(trg_file_path, 'w') as trg_file_opened:
                                                                if self.trg_file_type == 'rsids':
                                                                        trg_file_opened.write(ucsc_header_line + '\n')
                                                                        trg_file_opened.write('#rsID\n')
                                                                        trg_file_opened.write(query_snp_data.id + '\n')
                                                                elif self.trg_file_type == 'tsv':
                                                                        trg_file_opened.write(ucsc_header_line + '\n')
                                                                        trg_file_opened.write('#' + '\t'.join(header_row) + '\n')
                                                                        trg_file_opened.write('\t'.join(map(str, query_snp_ann)) + '\n')
                                                                elif self.trg_file_type == 'json':
                                                                        trg_obj = [dict(zip(meta_keys, meta_vals)),
                                                                                   dict(zip(header_row, query_snp_ann))]
                                                                        
                                                                #Перебор 1000-Genomes-SNPs
                                                                #в пределах окна, центром которого
                                                                #служит запрашиваемый SNP.
                                                                for oppos_snp_data in intgen_vcf_opened.fetch(chrom,
                                                                                                              low_bound,
                                                                                                              high_bound):
                                                                        
                                                                        #Оппонирующий SNP будет
                                                                        #отсеян, если совпадает с
                                                                        #запрашиваемым или имеет
                                                                        #не-rs-идентификатор или
                                                                        #является мультиаллельным.
                                                                        if oppos_snp_data.id == query_snp_data.id \
                                                                           or re.match(r'rs\d+$', oppos_snp_data.id) is None \
                                                                           or 'MULTI_ALLELIC' in oppos_snp_data.info:
                                                                                continue
                                                                        
                                                                        #Для SNP, прошедшего эту
                                                                        #фильтрацию, произведём
                                                                        #поиск генотипов по сэмплам.
                                                                        oppos_snp_genotypes = []
                                                                        for sample_name in self.sample_names:
                                                                                try:
                                                                                        oppos_snp_genotypes += oppos_snp_data.samples[sample_name]['GT']
                                                                                except KeyError:
                                                                                        continue
                                                                                
                                                                        #Получение значений LD с
                                                                        #помощью оффлайн-калькулятора.
                                                                        #Полезный побочный продукт функции -
                                                                        #частоты альтернативного аллеля
                                                                        #запрашиваемого и оппонирующего SNPs.
                                                                        trg_vals = calc_ld(query_snp_genotypes,
                                                                                           oppos_snp_genotypes)
                                                                        
                                                                        #Ещё один этап отбора оппонирующих
                                                                        #SNP будет по устновленному
                                                                        #исследователем порогу LD.
                                                                        if trg_vals[self.ld_thres_measure] < self.ld_low_thres:
                                                                                continue
                                                                        
                                                                        #Теперь понятно, что результаты
                                                                        #будут, значит, есть смысл
                                                                        #дать сигнал, запрещающий
                                                                        #удаление конечного файла.
                                                                        empty_res = False
                                                                        
                                                                        #Добавляем в конечный файл или объект
                                                                        #информацию о текущем оппонирующем SNP.
                                                                        #Если задан минималистический вывод,
                                                                        #то это будет только rsID, а если TSV
                                                                        #или JSON, то rsID с характеристиками.
                                                                        if self.trg_file_type == 'rsids':
                                                                                trg_file_opened.write(oppos_snp_data.id + '\n')
                                                                                continue
                                                                        oppos_snp_ann = [oppos_snp_data.pos,
                                                                                         oppos_snp_data.id,
                                                                                         oppos_snp_data.ref,
                                                                                         ','.join(oppos_snp_data.alts),
                                                                                         ','.join(oppos_snp_data.info['VT']),
                                                                                         trg_vals['snp_2_alt_freq'],
                                                                                         trg_vals['r_square'],
                                                                                         trg_vals['d_prime'],
                                                                                         oppos_snp_data.pos - query_snp_data.pos]
                                                                        if self.trg_file_type == 'tsv':
                                                                                trg_file_opened.write('\t'.join(map(str, oppos_snp_ann)) + '\n')
                                                                        elif self.trg_file_type == 'json':
                                                                                trg_obj.append(dict(zip(header_row, oppos_snp_ann)))
                                                                                
                                                                #При подготовке JSON-вывода
                                                                #результаты должны были накапливаться
                                                                #в объект - список словарей.
                                                                #Теперь надо прописать его в файл.
                                                                if self.trg_file_type == 'json':
                                                                        json.dump(trg_obj, trg_file_opened, indent=4)
                                                                        
                                                        #Если флаг-индикатор так и
                                                        #остался равен True, значит,
                                                        #результатов нет, и в конечный
                                                        #файл попали только хэдеры.
                                                        #Такие конечные файлы
                                                        #программа удалит.
                                                        if empty_res:
                                                                os.remove(trg_file_path)
                                                                
####################################################################################################

import sys, os, re, sqlite3, copy, json

#Подавление формирования питоновского кэша с
#целью предотвращения искажения результатов.
sys.dont_write_bytecode = True

from argparse import ArgumentParser, RawTextHelpFormatter
from backend.prep_intgen_data import prep_intgen_data
from backend.get_sample_names import get_sample_names
from multiprocessing import Pool
from pysam import VariantFile
from backend.calc_ld import calc_ld

#Подготовительный этап: обработка
#аргументов командной строки,
#создание экземпляра содержащего
#ключевую функцию класса, определение
#оптимального количества процессов.
args = add_args()
prep_single_proc = PrepSingleProc(args)
max_proc_quan = args.max_proc_quan
src_file_names = os.listdir(prep_single_proc.src_dir_path)
src_files_quan = len(src_file_names)
if max_proc_quan > src_files_quan <= 8:
        proc_quan = src_files_quan
elif max_proc_quan > 8:
        proc_quan = 8
else:
        proc_quan = max_proc_quan
        
print(f'\nПоиск SNPs в неравновесии по сцеплению')
print(f'\tколичество параллельных процессов: {proc_quan}')

#Параллельный запуск создания коллекций.
with Pool(proc_quan) as pool_obj:
        pool_obj.map(prep_single_proc.get_lnkd_snps, src_file_names)
