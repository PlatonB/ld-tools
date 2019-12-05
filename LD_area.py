__version__ = 'V7.2'

print('''
Программа ищет в пределах фланков SNPs,
обладающие надпороговым значением неравновесия
по сцеплению с каждым запрашиваемым SNP.

Автор: Платон Быкадоров (platon.work@gmail.com), 2018-2019.
Версия: V7.2.
Лицензия: GNU General Public License version 3.
Поддержать проект: https://money.yandex.ru/to/41001832285976
Документация: https://github.com/PlatonB/ld-tools/blob/master/README.md

Обязательно!
Перед запуском программы нужно установить модули:
pip3 install plyvel pysam --user

Поддерживаемые исходные файлы - таблицы,
содержащие столбец с набором refSNPIDs.
Если таких столбцов - несколько,
программа будет использовать самый левый.

Если настройки, запрашиваемые в рамках интерактивного
диалога, вам непонятны - пишите, пожалуйста, в Issues.
''')

def join_header_element(header_key, header_val):
        '''
        Конечные не-JSON-файлы будут начинаться
        с хэдеров, структура которых напоминает
        таковую у хэдеров UCSC Table Browser.
        Эта функция собирает каждый элемент хэдера.
        '''
        if type(header_val).__name__ == 'str':
                header_val = f'"{header_val}"'
        elif type(header_val).__name__ == 'list':
                header_val = ','.join([f'"{header_val_element}"' for header_val_element in header_val])
        return f'{header_key}={header_val}'

####################################################################################################

print('\nИмпорт модулей программы...')

import sys

#Подавление формирования питоновского кэша с
#целью предотвращения искажения результатов.
sys.dont_write_bytecode = True

import os, re, gzip, plyvel, json
from backend.prepare_intgen_data import process_intgen_data
from backend.retrieve_sample_indices import retrieve_sample_indices
from backend.ld_calc import ld_calc
from pysam import VariantFile

src_dir_path = input('\nПуть к папке с исходными файлами: ')

trg_top_dir_path = input('\nПуть к папке для конечных файлов: ')

intgen_dir_path = input('\nПуть к папке для данных 1000 Genomes: ')

num_of_headers = input('''\nКоличество не обрабатываемых строк
в начале каждой исходной таблицы
(игнорирование ввода ==> хэдеров/шапок в таблицах нет)
[0(|<enter>)|1|2|...]: ''')
if num_of_headers == '':
        num_of_headers = 0
else:
        num_of_headers = int(num_of_headers)
        
populations = input('''\nДля индивидов какой(-их) супер-/субпопуляции(-ий) считать LD?
(http://www.internationalgenome.org/faq/which-populations-are-part-your-study/)
(несколько - через пробел)
(игнорирование ввода ==> всех)
[all(|<enter>)|eur|jpt|jpt amr yri|...]: ''').upper().split()
if populations == []:
        populations = ['ALL']
        
genders = input('''\nДля индивидов какого пола считать LD?
(игнорирование ввода ==> обоих)
[male(|m)|female(|f)|both(|<enter>)]: ''').split()
if genders == ['m']:
        genders = ['male']
elif genders == ['f']:
        genders = ['female']
elif genders in [[], ['both']]:
        genders = ['male', 'female']
elif genders not in [['male'], ['female']]:
        print(f'{" ".join(genders)} - недопустимая опция')
        sys.exit()
        
flank_size = input('''\nРазмер *каждого* из фланков, в пределах
которых надо искать сцепленные SNP
(игнорирование ввода ==> 100000)
[100000(|<enter>)|250000|...]: ''')
if flank_size == '':
        flank_size = 100000
else:
        flank_size = int(flank_size)
        
thres_ld_measure = input('''\nМера LD для выставления порога
[r_square(|r)|d_prime(|d)]: ''')
if thres_ld_measure == 'r':
        thres_ld_measure = 'r_square'
elif thres_ld_measure == 'd':
        thres_ld_measure = 'd_prime'
elif thres_ld_measure not in ['r_square', 'd_prime']:
        print(f'{thres_ld_measure} - недопустимая опция')
        sys.exit()
        
thres_ld_value = float(input(f'\n{thres_ld_measure} ≥ '))

trg_file_type = input('''\nФормат конечных файлов
(игнорирование ввода ==> json)
[json(|j|<enter>)|tsv(|t)|rsids(|r)]: ''')
if trg_file_type in ['j', '']:
        trg_file_type = 'json'
elif trg_file_type == 't':
        trg_file_type = 'tsv'
elif trg_file_type == 'r':
        trg_file_type = 'rsids'
elif trg_file_type not in ['json', 'tsv', 'rsids']:
        print(f'{trg_file_type} - недопустимая опция')
        sys.exit()
        
verbose = input('''\nВыводить подробную информацию о ходе работы?
(не рекомендуется, если набор
исходных данных очень большой)
(игнорирование ввода ==> не выводить)
[yes(|y)|no(|n|<enter>)]: ''')
if verbose not in ['yes', 'y', 'no', 'n', '']:
        print(f'{verbose} - недопустимая опция')
        sys.exit()
        
#Вызов функции, которая, во-первых, скачает
#заархивированные VCF и панель сэмплов проекта
#1000 Genomes, если они ещё не скачаны, во-вторых,
#сконвертирует панель в JSON, а данные изо
#всех VCF разместит в единую LevelDB-базу.
#LevelDB-база данных (далее - база) будет
#состоять из ключей - ID снипов всех хромосом,
#и значений - сжатых строк упомянутых VCF.
intgen_sampjson_path, intgen_vcfdb_path, intgen_vcfgz_paths = process_intgen_data(intgen_dir_path)

#Вызов функции, производящей отбор
#сэмплов, относящихся к указанным
#пользователем популяциям и полам,
#и возвращающей индексы этих сэмплов.
sample_indices = retrieve_sample_indices(intgen_sampjson_path,
                                         populations,
                                         genders,
                                         intgen_vcfdb_path)

#Открытие базы.
intgen_vcfdb_opened = plyvel.DB(intgen_vcfdb_path, create_if_missing=False)

#Работа с исходными файлами.
src_file_names = os.listdir(src_dir_path)
for src_file_name in src_file_names:
        if src_file_name.startswith('.~lock.'):
                continue
        src_file_base = '.'.join(src_file_name.split('.')[:-1])
        
        #Построение пути к папке для результатов
        #по текущему пользовательскому файлу,
        #которые, в свою очередь, должны будут
        #распределяться по хромосомным подпапкам.
        trg_dir_path = os.path.join(trg_top_dir_path, f'{src_file_base}_lnkd')
        
        #Открытие файла пользователя на чтение.
        with open(os.path.join(src_dir_path, src_file_name)) as src_file_opened:
                
                print(f'''\n~
\n{src_file_name}''')
                
                #Считываем строки-хэдеры, чтобы сместить
                #курсор к началу основной части таблицы.
                for header_index in range(num_of_headers):
                        src_file_opened.readline()
                        
                #Построчное прочтение основной части исходной таблицы.
                for line in src_file_opened:
                        
                        #Попытка получения refSNPID из текущей строки.
                        try:
                                query_rs_id = re.search(r'rs\d+', line).group()
                        except AttributeError:
                                continue
                        
                        if verbose in ['yes', 'y']:
                                print(f'\n{query_rs_id}...')
                                
                        #Список, в который будут помещаться,
                        #как минимум, идентификаторы SNP,
                        #сцепленных с запрашиваемым,
                        #и объект-предшественник
                        #стартового элемента (словаря
                        #или хэдера с метаинформацией).
                        #Туда также могут попадать
                        #различные характеристики SNPs.
                        #Содержимое списка будет зависить
                        #от пресета, выбранного пользователем.
                        linked_snps = []
                        
                        #Попытка извлечения из очередной базы сжатой строки,
                        #соответствующей текущему пользовательскому refSNPID.
                        #Декомпрессия результата, превращение из
                        #байт-строки в обычную и преобразование в список.
                        #Если refSNPID в базе нет, то вместо сжатой строки вернётся None.
                        #В результате разархивации None выйдет пустая байтовая строка.
                        #После её конвертаций, упомянутых выше, сформируется список
                        #с обычной пустой строкой в качестве единственного элемента.
                        query_snp_row = gzip.decompress(intgen_vcfdb_opened.get(query_rs_id.encode())).decode('utf-8').split('\t')
                        
                        #Если refSNPID из файла
                        #пользователя не обнаружился
                        #в базе, значит он - невалидный.
                        #В матрицу он допущен не будет.
                        #Одна из вероятных причин - в том,
                        #что именуемый им SNP - не биаллельный.
                        #Переходим к следующей строке.
                        if query_snp_row == ['']:
                                if verbose in ['yes', 'y']:
                                        print('\tневалидный refSNPID (возможно, ID мультиаллельного SNP).')
                                continue
                                
                        #Получение номера хромосомы запрашиваемого
                        #и впоследствии находимых SNP.
                        chr_num, query_snp_pos = query_snp_row[0], \
                                                 int(query_snp_row[1])
                        
                        #Путь к хромосомной подпапке, имя
                        #конечного файла и путь к нему.
                        #Не факт, что соответствующие подпапка
                        #и файл будут далее реально созданы.
                        #Пути не пригодятся, если упомянутые
                        #объекты файловой системы появились в
                        #одной из прошлых итераций цикла.
                        #Если данных объектов ещё нет, но
                        #для текущего SNP не нашлось ни одного
                        #сцепленного, то в файловой системе
                        #также никаких изменений не произойдёт.
                        trg_chrdir_path = os.path.join(trg_dir_path, chr_num)
                        if trg_file_type in ['json', 'tsv']:
                                ext = trg_file_type
                        else:
                                ext = 'txt'
                        trg_file_name = f'{chr_num}_{query_rs_id}_{thres_ld_measure[0]}_{str(thres_ld_value)}.{ext}'
                        trg_file_path = os.path.join(trg_chrdir_path, trg_file_name)
                        
                        #Если запрашиваемый SNP встретился
                        #повторно, то необходимо предотвратить
                        #поиск SNPs, с ним сцепленных.
                        if os.path.exists(trg_file_path) == True:
                                if verbose in ['yes', 'y']:
                                        print('\tуже был ранее обработан')
                                continue
                        
                        #Перебор путей к 1000 Genomes-архивам и поиск
                        #того архива, в котором спрятан запрашиваемый SNP.
                        #Поиск остановится после обнаружения вхождения
                        #шаблона, содержащего ранее полученный номер
                        #хромосомы запрашиваемого SNP, в имя архива.
                        for intgen_vcfgz_path in intgen_vcfgz_paths:
                                if os.path.basename(intgen_vcfgz_path).find(f'chr{chr_num}_') != -1:
                                        
                                        #Открытие на чтение того из tabix-индексированных
                                        #1000 Genomes-VCF, в котором точно есть запрашиваемый SNP.
                                        #Pysam требует на вход лишь путь к vcf.gz-файлу,
                                        #а tbi-сателлита обнаруживает самостоятельно.
                                        #Из этого vcf.gz-файла будут вытаскиваться сцепленные SNPs.
                                        with VariantFile(intgen_vcfgz_path) as variant_file_obj:
                                                
                                                #Получение координат — границ фланков, размер
                                                #которых на этапе диалога задан пользователем.
                                                #Расстояние между границами будем называть окном.
                                                #Если запрашиваемый SNP расположен настолько близко к
                                                #началу хромосомы, что левый фланк вылезает за него, то
                                                #окно будет принудительно отмеряться от нулевой координаты.
                                                #Аналогичной проблемы для конца хромосомы не возникает,
                                                #т.к. здесь pysam — молодец — самостоятельно справляется.
                                                window_left_border = query_snp_pos - flank_size
                                                if window_left_border < 0:
                                                        window_left_border = 0
                                                window_right_border = query_snp_pos + flank_size
                                                
                                                #Перебор объектов, соответствующих окну.
                                                #Каждый объект - это представленная
                                                #pysam'ом refSNPID-содержащая строка.
                                                #SNPs из описываемых объектов
                                                #будем называть кандидатными.
                                                for rec in variant_file_obj.fetch(chr_num,
                                                                                  window_left_border,
                                                                                  window_right_border):
                                                        
                                                        #Преобразование pysam-объекта
                                                        #в строку, а её - в список.
                                                        oppos_snp_row = str(rec).split('\n')[0].split('\t')
                                                        
                                                        #Получение из этого списка ID
                                                        #и info-ячейки кандидатного SNP.
                                                        #В последней обязательно присутствует
                                                        #уточнённый тип мутации, а также может
                                                        #встречаться флаг MULTI_ALLELIC, позволяющий
                                                        #идентифицировать не-биаллельные SNP.
                                                        oppos_rs_id, oppos_snp_info = oppos_snp_row[2], \
                                                                                      oppos_snp_row[7]
                                                        
                                                        #Кандидатный SNP имеет шанс войти в конечный
                                                        #список при соответствии нескольким критериям:
                                                        #его идентификатор - refSNPID, он - не мультиаллельный
                                                        #и не тот же самый, что и запрашиваемый.
                                                        if re.match(r'rs\d+$', oppos_rs_id) != None and \
                                                           oppos_snp_info.find('MULTI_ALLELIC') == -1 and \
                                                           oppos_rs_id != query_rs_id:
                                                                
                                                                #Получение значений LD.
                                                                #Полезный побочный продукт функции -
                                                                #частоты альтернативного аллеля
                                                                #запрашиваемого и кандидатного SNPs.
                                                                trg_vals = ld_calc(query_snp_row,
                                                                                   oppos_snp_row,
                                                                                   sample_indices)
                                                                
                                                                #Отбор кандидатных SNP по
                                                                #пользовательскому порогу LD.
                                                                if trg_vals[thres_ld_measure] < thres_ld_value:
                                                                        continue
                                                                
                                                                #Добавление в конечный список очередного элемента.
                                                                #Что это будет за элемент - зависит от решения
                                                                #пользователя, в каком виде выводить результаты.
                                                                #Для начала проверим, не предпочёл ли пользователь
                                                                #самый минималистичный вывод - набор refSNPID.
                                                                #Если да, заселим в список текущий refSNPID и сразу
                                                                #перейдём к строке со следующим кандидатным SNP.
                                                                if trg_file_type == 'rsids':
                                                                        linked_snps.append(oppos_rs_id)
                                                                        continue
                                                                
                                                                #Если конечный формат - JSON
                                                                #или TSV, то тут всё сложнее.
                                                                #Файлы таких форматов для того
                                                                #и создаём, чтобы в них уместить
                                                                #в структурированном виде
                                                                #как можно больше информации.
                                                                #Ради дальнейшего размещения
                                                                #в конечный список соберём
                                                                #позицию, ID, оба сиквенса,
                                                                #точный тип мутации и частоту
                                                                #альтернативного аллеля
                                                                #(AF) текущего сцепленного
                                                                #SNP, а также значения LD и
                                                                #физическое расстояние между
                                                                #этим SNP и запрашиваемым.
                                                                oppos_snp_pos, oppos_snp_ref, oppos_snp_alt = int(oppos_snp_row[1]), \
                                                                                                              oppos_snp_row[3], \
                                                                                                              oppos_snp_row[4]
                                                                oppos_snp_vals = [oppos_snp_pos,
                                                                                  oppos_rs_id,
                                                                                  oppos_snp_ref,
                                                                                  oppos_snp_alt,
                                                                                  re.search(r'(?<=VT=)\w+?\b', oppos_snp_info).group(),
                                                                                  trg_vals['snp_2_alt_freq'],
                                                                                  trg_vals['r_square'],
                                                                                  trg_vals['d_prime'],
                                                                                  oppos_snp_pos - query_snp_pos]
                                                                
                                                                #Если запланирован JSON-вывод,
                                                                #для собранных значений составим
                                                                #список соответствующих ключей.
                                                                #Объединим ключи и значения в словарь
                                                                #и добавим его в конечный список.
                                                                if trg_file_type == 'json':
                                                                        oppos_snp_keys = ['lnkd.hg38_pos',
                                                                                          'lnkd.rsID',
                                                                                          'lnkd.ref',
                                                                                          'lnkd.alt',
                                                                                          'lnkd.type',
                                                                                          'lnkd.alt_freq',
                                                                                          "r2",
                                                                                          "D'",
                                                                                          'dist']
                                                                        linked_snps.append(dict(zip(oppos_snp_keys, oppos_snp_vals)))
                                                                        
                                                                #Если output - TSV - упомянутые значения
                                                                #сконкатенируются в табличную строку,
                                                                #которая поступит в конечный список.
                                                                elif trg_file_type == 'tsv':
                                                                        linked_snps.append('\t'.join([str(oppos_snp_val) for oppos_snp_val in oppos_snp_vals]))
                                                                        
                                        #Конечная папка, хромосомная подпапка и
                                        #файл с результатами могут быть созданы
                                        #при условии, что в конечном списке
                                        #оказался хоть один сцепленный SNP.
                                        if linked_snps != []:
                                                
                                                #Создание конечной папки и хромосомной
                                                #подпапки, если таковых ещё нет.
                                                #Если подпапка есть, то наличие
                                                #папки проверяться не будет.
                                                if os.path.exists(trg_chrdir_path) == False:
                                                        if os.path.exists(trg_dir_path) == False:
                                                                os.mkdir(trg_dir_path)
                                                        os.mkdir(trg_chrdir_path)
                                                        
                                                #Если файлу с результатами быть, подготовим
                                                #список таких метазначений, как номер
                                                #хромосомы и параметры запроса, плюс список
                                                #названий каждого из них (далее - метаключей).
                                                #Метазначения с метаключами в том или ином виде
                                                #должны будут пойти в конечный файл любого формата.
                                                common_ann_keys, common_ann_vals = ['chr',
                                                                                    'each_flank',
                                                                                    f'{thres_ld_measure}_thres',
                                                                                    'pops',
                                                                                    'gends'], [chr_num,
                                                                                               flank_size,
                                                                                               thres_ld_value,
                                                                                               populations,
                                                                                               genders]
                                                
                                                #Для дальнейшего размещения в JSON или TSV
                                                #будут также собраны сведения о запрашиваемом
                                                #SNP: позиция, ID, сиквенсы, конкретизированный
                                                #тип, плюс частота альтернативного аллеля.
                                                #AF запрашиваемого SNP возвращается
                                                #при каждом вызове функции расчёта
                                                #LD, поэтому её можно взять из словаря,
                                                #получившегося при последнем таком расчёте.
                                                #Ещё что нужно для JSON и таблицы, так это
                                                #ключи, представляющие собой названия
                                                #характеристик того или иного SNP.
                                                if trg_file_type in ['json', 'tsv']:
                                                        query_snp_ref, query_snp_alt, query_snp_type = query_snp_row[3], \
                                                                                                       query_snp_row[4], \
                                                                                                       re.search(r'(?<=VT=)\w+?\b', query_snp_row[7]).group()
                                                        query_ann_vals, header_keys = [query_snp_pos,
                                                                                       query_rs_id,
                                                                                       query_snp_ref,
                                                                                       query_snp_alt,
                                                                                       query_snp_type,
                                                                                       trg_vals['snp_1_alt_freq']], ['hg38_pos',
                                                                                                                     'rsID',
                                                                                                                     'ref',
                                                                                                                     'alt',
                                                                                                                     'type',
                                                                                                                     'alt_freq']
                                                        
                                                #В конечный список, собираемый
                                                #под таблицу или refSNPIDs-файл,
                                                #добавим хэдер сантакрузовского
                                                #стиля, построив его из ранее
                                                #описанных метаключей и метазначений.
                                                if trg_file_type in ['tsv', 'rsids']:
                                                        linked_snps.insert(0, '##' + ' '.join(map(join_header_element,
                                                                                                 common_ann_keys,
                                                                                                 common_ann_vals)))
                                                        
                                                #Создание и открытие конечного файла на запись.
                                                with open(trg_file_path, 'w') as trg_file_opened:
                                                        
                                                        #JSON: первый словарь получаем из уже
                                                        #знакомых нам метаключей и метазначений,
                                                        #второй - из ключей, к каждому из которых
                                                        #добавляется приставка, указывающая на
                                                        #призвание пары описывать запрашиваемый
                                                        #SNP, и соответствующих значений.
                                                        #Размещаем список словарей в файл,
                                                        #автоматически проставляя отступы.
                                                        #Телепортируемся через огромную толщу
                                                        #кода к следующей строке исходного файла:).
                                                        if trg_file_type == 'json':
                                                                linked_snps.insert(0, dict(zip(common_ann_keys,
                                                                                               common_ann_vals)))
                                                                linked_snps.insert(1, dict(zip([f'quer.{header_key}' for header_key in header_keys],
                                                                                               query_ann_vals)))
                                                                json.dump(linked_snps, trg_file_opened, indent=4)
                                                                break
                                                        
                                                        #TSV: Собираем, во-первых, табулированную шапку
                                                        #из названий сниповых характеристик, добавив
                                                        #к ним элементы r2, D' и dist, а, во-вторых,
                                                        #строку с характеристиками запрашиваемого SNP.
                                                        elif trg_file_type == 'tsv':
                                                                linked_snps.insert(1, '#' + '\t'.join(header_keys + ["r2",
                                                                                                                     "D'",
                                                                                                                     'dist']))
                                                                linked_snps.insert(2, '\t'.join([str(query_ann_val) for query_ann_val in query_ann_vals] + ['quer'] * 3))
                                                                
                                                        #RefSNPIDs: дополняем конечный
                                                        #список именем снипового столбца и
                                                        #идентификатором запрашиваемого SNP.
                                                        elif trg_file_type == 'rsids':
                                                                linked_snps.insert(1, '#rsID')
                                                                linked_snps.insert(2, query_rs_id)
                                                                
                                                        #Независимо от того, выбран табличный
                                                        #формат или вывод одиночных refSNPID,
                                                        #прописываем все собранные в конечный
                                                        #список данные как обычные строки.
                                                        for line in linked_snps:
                                                                trg_file_opened.write(line + '\n')
                                                                
                                        elif verbose in ['yes', 'y']:
                                                print(f'\t{thres_ld_measure} >= {thres_ld_value}: SNPs в LD не найдены')
                                                
                                        break
                                
intgen_vcfdb_opened.close()
