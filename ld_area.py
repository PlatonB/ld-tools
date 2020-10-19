__version__ = 'V11.1'

def build_ucsc_header(header_key, header_val):
        '''
        Конечные не-JSON-файлы будут начинаться
        с хэдеров, структура которых напоминает
        таковую у хэдеров UCSC Table Browser. Эта
        функция собирает каждый элемент хэдера.
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
        поиск вариантов, неравновесно
        сцепленных с вариантами
        из исходного файла.
        '''
        def __init__(self, args):
                '''
                Получение атрибутов, необходимых заточенной
                под многопроцессовое выполнение функции
                извлечения вариантов, обладающих надпороговым
                LD по отношению к вариантам исходного
                набора. Атрибуты должны быть созданы
                единожды и далее ни в коем случае
                не изменяться. Получаются они в основном
                из указанных исследователем опций.
                '''
                self.src_dir_path = os.path.normpath(args.src_dir_path)
                self.intgen_dir_path = os.path.normpath(args.intgen_dir_path)
                if args.trg_top_dir_path is None:
                        self.trg_top_dir_path = self.src_dir_path
                else:
                        self.trg_top_dir_path = os.path.normpath(args.trg_top_dir_path)
                self.meta_lines_quan = args.meta_lines_quan
                if args.skip_intgen_data_ver:
                        self.intgen_convdb_path = os.path.join(self.intgen_dir_path,
                                                               'conversion.db')
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
                
        def get_inld_vars(self, src_file_name):
                '''
                Функция нахождения вариантов, находящихся
                в надпороговом неравновесии по сцеплению
                с каждым вариантом исходной таблицы.
                '''
                
                #Считывание исходной таблицы, извлечение оттуда
                #rsIDs и создание словаря, в котором позиции и
                #идентификаторы вариантов разбиты по хромосомам.
                data_by_chrs = create_src_dict(self.src_dir_path,
                                               src_file_name,
                                               self.meta_lines_quan,
                                               self.intgen_convdb_path)
                
                #В одну папку второго уровня планируется
                #размещать все результаты, полученные по
                #данным одного исходного файла. Папка реально
                #создастся только тогда, когда появится
                #гарантия, что она не останется пустой.
                src_file_base = src_file_name.rsplit('.', maxsplit=1)[0]
                trg_dir_path = os.path.join(self.trg_top_dir_path,
                                            f'{src_file_base}_in_LD')
                
                #Расширение всех конечных файлов зависит
                #от выбранного исследователем формата.
                if self.trg_file_type in ['tsv', 'json']:
                        ext = self.trg_file_type
                else:
                        ext = 'txt'
                        
                #Сокращённые заглавия таких
                #глобальных характеристик,
                #как хромосома и настройки
                #отбора оппонирующих вариантов.
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
                for chrom in data_by_chrs:
                        chr_dir_path = os.path.join(trg_dir_path,
                                                    chrom)
                        os.makedirs(chr_dir_path)
                        
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
                        
                        #Поскольку известно, какая хромосома сейчас
                        #обрабатывается, можно сразу открыть архив
                        #именно с относящимися к этой хромосоме данными
                        #1000 Genomes. Из него для начала будем извлекать
                        #аннотации только запрашиваемых вариантов.
                        #Нужный ли вариант нашёлся в 1000 Genomes,
                        #дополнительно проверяется по rs-идентификатору.
                        with VariantFile(os.path.join(self.intgen_dir_path,
                                                      f'{chrom}.vcf.gz')) as intgen_vcf_opened:
                                for var_row in data_by_chrs[chrom]:
                                        for intgen_rec in intgen_vcf_opened.fetch(chrom,
                                                                                  var_row[0] - 1,
                                                                                  var_row[0]):
                                                if intgen_rec.id != var_row[1]:
                                                        continue
                                                query_var_rec = intgen_rec
                                                break
                                        trg_file_name = f'{query_var_rec.id}_chr{chrom}_{self.ld_thres_measure[0]}_{str(self.ld_low_thres)}.{ext}'
                                        trg_file_path = os.path.join(chr_dir_path, trg_file_name)
                                        
                                        #Создаём флаг, по которому далее будет
                                        #определено, оказались ли в конечном
                                        #файле строки, отличные от хэдеров.
                                        empty_res = True
                                        
                                        #Вычисляем границы фланков вокруг текущего
                                        #запрашиваемого варианта, где надо искать
                                        #неравновесно сцепленные с ним варианты.
                                        #Если нижняя граница окна заползёт
                                        #за ноль, то во избежание ошибки со
                                        #стороны pysam приравняем её к нулю.
                                        low_bound = query_var_rec.pos - self.flank_size
                                        if low_bound < 0:
                                                low_bound = 0
                                        high_bound = query_var_rec.pos + self.flank_size
                                        
                                        #Получение генотипов и
                                        #ключевых характеристик
                                        #запрашиваемого варианта.
                                        query_var_genotypes = []
                                        for sample_name in self.sample_names:
                                                try:
                                                        query_var_genotypes += query_var_rec.samples[sample_name]['GT']
                                                except KeyError:
                                                        continue
                                        query_var_alt_freq = round(query_var_genotypes.count(1) /
                                                                   len(query_var_genotypes), 4)
                                        query_var_ann = [query_var_rec.pos,
                                                         query_var_rec.id,
                                                         query_var_rec.ref,
                                                         ','.join(query_var_rec.alts),
                                                         ','.join(query_var_rec.info['VT']),
                                                         query_var_alt_freq] + ['quer'] * 3
                                        
                                        #Открываем конечный файл на запись
                                        #и проделываем кропотливую работу по
                                        #созданию формат-ориентированных хэдеров.
                                        with open(trg_file_path, 'w') as trg_file_opened:
                                                if self.trg_file_type == 'rsids':
                                                        trg_file_opened.write(ucsc_header_line + '\n')
                                                        trg_file_opened.write('#rsID\n')
                                                        trg_file_opened.write(query_var_rec.id + '\n')
                                                elif self.trg_file_type == 'tsv':
                                                        trg_file_opened.write(ucsc_header_line + '\n')
                                                        trg_file_opened.write('#' + '\t'.join(header_row) + '\n')
                                                        trg_file_opened.write('\t'.join(map(str, query_var_ann)) + '\n')
                                                elif self.trg_file_type == 'json':
                                                        trg_obj = [dict(zip(meta_keys, meta_vals)),
                                                                   dict(zip(header_row, query_var_ann))]
                                                        
                                                #Перебор 1000-Genomes-вариантов в пределах окна,
                                                #центром которого служит запрашиваемый вариант.
                                                for oppos_var_rec in intgen_vcf_opened.fetch(chrom,
                                                                                             low_bound,
                                                                                             high_bound):
                                                        
                                                        #Оппонирующий вариант будет отсеян, если
                                                        #совпадает с запрашиваемым или имеет не-rs-
                                                        #идентификатор или является мультиаллельным.
                                                        if oppos_var_rec.id == query_var_rec.id \
                                                           or re.match(r'rs\d+$', oppos_var_rec.id) is None \
                                                           or 'MULTI_ALLELIC' in oppos_var_rec.info:
                                                                continue
                                                        
                                                        #Для варианта, прошедшего
                                                        #эту фильтрацию, произведём
                                                        #поиск генотипов по сэмплам.
                                                        oppos_var_genotypes = []
                                                        for sample_name in self.sample_names:
                                                                try:
                                                                        oppos_var_genotypes += oppos_var_rec.samples[sample_name]['GT']
                                                                except KeyError:
                                                                        continue
                                                                
                                                        #Получение значений LD с помощью
                                                        #оффлайн-калькулятора. Полезный
                                                        #побочный продукт функции - частоты
                                                        #альтернативного аллеля запрашиваемого
                                                        #и оппонирующего вариантов.
                                                        trg_vals = calc_ld(query_var_genotypes,
                                                                           oppos_var_genotypes)
                                                        
                                                        #Ещё один этап отбора оппонирующих
                                                        #вариантов будет по установленному
                                                        #исследователем порогу LD.
                                                        if trg_vals[self.ld_thres_measure] < self.ld_low_thres:
                                                                continue
                                                        
                                                        #Теперь понятно, что результаты
                                                        #будут, значит, есть смысл
                                                        #дать сигнал, запрещающий
                                                        #удаление конечного файла.
                                                        empty_res = False
                                                        
                                                        #Добавляем в конечный файл или объект информацию
                                                        #о текущем оппонирующем варианте. Если задан
                                                        #минималистический вывод, то это будет только rsID,
                                                        #а если TSV или JSON, то rsID с характеристиками.
                                                        if self.trg_file_type == 'rsids':
                                                                trg_file_opened.write(oppos_var_rec.id + '\n')
                                                                continue
                                                        oppos_var_ann = [oppos_var_rec.pos,
                                                                         oppos_var_rec.id,
                                                                         oppos_var_rec.ref,
                                                                         ','.join(oppos_var_rec.alts),
                                                                         ','.join(oppos_var_rec.info['VT']),
                                                                         trg_vals['var_2_alt_freq'],
                                                                         trg_vals['r_square'],
                                                                         trg_vals['d_prime'],
                                                                         oppos_var_rec.pos - query_var_rec.pos]
                                                        if self.trg_file_type == 'tsv':
                                                                trg_file_opened.write('\t'.join(map(str, oppos_var_ann)) + '\n')
                                                        elif self.trg_file_type == 'json':
                                                                trg_obj.append(dict(zip(header_row, oppos_var_ann)))
                                                                
                                                #При подготовке JSON-вывода
                                                #результаты должны были накапливаться
                                                #в список словарей. Теперь
                                                #надо прописать его в файл.
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

import sys, locale, os, re, json

#Подавление формирования питоновского кэша с
#целью предотвращения искажения результатов.
sys.dont_write_bytecode = True

from cli.ld_area_cli_ru import add_args_ru
from cli.ld_area_cli_en import add_args_en
from backend.prep_intgen_data import prep_intgen_data
from backend.get_sample_names import get_sample_names
from backend.create_src_dict import create_src_dict
from multiprocessing import Pool
from pysam import VariantFile
from backend.calc_ld import calc_ld

#Подготовительный этап: обработка
#аргументов командной строки,
#создание экземпляра содержащего
#ключевую функцию класса, определение
#оптимального количества процессов.
if locale.getdefaultlocale()[0][:2] == 'ru':
        args = add_args_ru(__version__)
else:
        args = add_args_en(__version__)
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
        
print(f'\nSelecting variants in LD and in window')
print(f'\tnumber of parallel processes: {proc_quan}')

#Параллельный запуск создания коллекций.
with Pool(proc_quan) as pool_obj:
        pool_obj.map(prep_single_proc.get_inld_vars,
                     src_file_names)
