__version__ = 'V8.3'

print('''
Программа ищет в пределах фланков SNPs,
обладающие надпороговым значением неравновесия
по сцеплению с каждым запрашиваемым SNP.

Автор: Платон Быкадоров (platon.work@gmail.com), 2018-2020.
Версия: V8.3.
Лицензия: GNU General Public License version 3.
Поддержать проект: https://money.yandex.ru/to/41001832285976
Документация: https://github.com/PlatonB/ld-tools/blob/master/README.md

Обязательно!
Перед запуском программы нужно установить
MongoDB и PyMongo (см. README).

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

import sys, os, re, json

#Подавление формирования питоновского кэша с
#целью предотвращения искажения результатов.
sys.dont_write_bytecode = True

from backend.init_intgen_db import *
from backend.create_intgen_db import create_intgen_db
from backend.samp_manager import get_sample_names, get_phased_genotypes
from backend.ld_calc import ld_calc

src_dir_path = input('\nПуть к папке с исходными файлами: ')

trg_top_dir_path = input('\nПуть к папке для конечных файлов: ')

num_of_headers = input('''\nКоличество не обрабатываемых строк
в начале каждой исходной таблицы
(игнорирование ввода ==> хэдеров/шапок в таблицах нет)
[0(|<enter>)|1|2|...]: ''')
if num_of_headers == '':
        num_of_headers = 0
else:
        num_of_headers = int(num_of_headers)
        
populations = input('''\nДля индивидов каких супер-/субпопуляций считать LD?
(http://www.internationalgenome.org/faq/which-populations-are-part-your-study/)
(несколько - через пробел)
[all(|<enter>)|eur|jpt|jpt amr yri|...]: ''').upper().split()
if populations == []:
        populations = ['ALL']
        
genders = input('''\nДля индивидов какого пола считать LD?
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
[yes(|y)|no(|n|<enter>)]: ''')
if verbose not in ['yes', 'y', 'no', 'n', '']:
        print(f'{verbose} - недопустимая опция')
        sys.exit()
        
#Отдельным модулем ранее были созданы
#объекты базы данных и коллекции сэмплов.
#Но существование объектов не означает, что
#база уже заполнена данными 1000 Genomes.
#Проверяем это по наличию созданных коллекций.
#Если нет ни одной коллекции, то запускаем
#модуль построения 1000 Genomes-базы с нуля.        
if intgen_db_obj.command('dbstats')['collections'] == 0:
        create_intgen_db()
        
#Вызов функции, производящей отбор
#сэмплов, относящихся к указанным
#исследователем популяциям и полам.
sample_names = get_sample_names(populations,
                                genders)

#Получение списка имён всех коллекций.
intgen_coll_names = intgen_db_obj.list_collection_names()

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
        
        #Открытие файла исследователя на чтение.
        with open(os.path.join(src_dir_path, src_file_name)) as src_file_opened:
                
                print(f'''\n~
\n{src_file_name}''')
                
                #Считываем строки-хэдеры, чтобы сместить
                #курсор к началу основной части таблицы.
                for header_index in range(num_of_headers):
                        src_file_opened.readline()
                        
                #Построчное прочтение основной
                #части исходной таблицы.
                for line in src_file_opened:
                        
                        #Попытка получения refSNPID из текущей строки.
                        try:
                                query_snp_id = re.search(r'rs\d+', line).group()
                        except AttributeError:
                                continue
                        
                        if verbose in ['yes', 'y']:
                                print(f'\n{query_snp_id}')
                                
                        #База была спроектирована так,
                        #чтобы в каждой коллекции
                        #были данные одной хромосомы.
                        #Коллекцию сэмплов я здесь,
                        #разумеется, не имею в виду.
                        #Текущий refSNPID из файла
                        #исследователя будет
                        #искаться поочерёдно
                        #в хромосомных коллекциях.
                        #Поиск будет прерван в
                        #случае обнаружения вхождения.
                        #Успешный результат - словарь с
                        #refSNPID, характеристиками SNP
                        #и фазированными генотипами.
                        #Неуспешный - None-значение.
                        #В первом случае мы сохраняем
                        #не только словарь, но и объект
                        #коллекции, в которой нашёлся
                        #текущий сниповый идентификатор.
                        #В последнем случае этот refSNPID
                        #объявляется невалидным, и
                        #осуществляется переход к
                        #обработке следующего refSNPID.
                        for intgen_coll_name in intgen_coll_names:
                                if intgen_coll_name == 'samples':
                                        continue
                                query_snp_data = intgen_db_obj[intgen_coll_name].find_one({'ID': query_snp_id})
                                if query_snp_data != None:
                                        intgen_vcfcoll_obj = intgen_db_obj[intgen_coll_name]
                                        query_snp_phased_genotypes = get_phased_genotypes(query_snp_data, sample_names)
                                        break
                        else:
                                if verbose in ['yes', 'y']:
                                        print('\tневалидный refSNPID (возможно, ID мультиаллельного SNP).')
                                continue
                                
                        #Получение номера хромосомы запрашиваемого
                        #и впоследствии находимых SNP.
                        chrom = query_snp_data['#CHROM']
                        
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
                        trg_chrdir_path = os.path.join(trg_dir_path, chrom)
                        if trg_file_type in ['json', 'tsv']:
                                ext = trg_file_type
                        else:
                                ext = 'txt'
                        trg_file_name = f'{chrom}_{query_snp_id}_{thres_ld_measure[0]}_{str(thres_ld_value)}.{ext}'
                        trg_file_path = os.path.join(trg_chrdir_path, trg_file_name)
                        
                        #Если запрашиваемый SNP встретился
                        #повторно, то необходимо предотвратить
                        #поиск SNPs, с ним сцепленных.
                        if os.path.exists(trg_file_path) == True:
                                if verbose in ['yes', 'y']:
                                        print('\tуже был ранее обработан')
                                continue
                        
                        #Вытаскиваем из найденного
                        #в коллекции словаря позицию
                        #текущего запрашиваемого SNP.
                        #Создаём список, в который
                        #будут помещаться, как минимум,
                        #идентификаторы SNP, сцепленных с
                        #запрашиваемым, и объект-предшественник
                        #стартового элемента (словаря
                        #или хэдера с метаинформацией).
                        #Туда также могут попадать
                        #различные характеристики SNPs.
                        #Содержимое списка будет зависить
                        #от пресета, выбранного исследователем.
                        query_snp_pos, linked_snps = query_snp_data['POS'], []
                        
                        #Позиция запрашиваемого SNP есть,
                        #теперь можно вытащить из той же
                        #коллекции, где был найден этот SNP,
                        #словари со всеми SNP в заданной
                        #исследователем окрестности.
                        #Далее эти словари и соответствующие
                        #SNP будут называться кандидатными.
                        oppos_curs_obj = intgen_vcfcoll_obj.find({'POS': {'$gte': query_snp_pos - flank_size,
                                                                          '$lte': query_snp_pos + flank_size}})
                        
                        #Перебор всех кандидатных словарей,
                        #оппонирующих текущему refSNPID.
                        #Среди кандидатных идентификаторов
                        #обязательно попадётся запрашиваемый.
                        #Он, разумеется, будет отсеян.
                        for oppos_snp_data in oppos_curs_obj:
                                oppos_snp_id = oppos_snp_data['ID']
                                if oppos_snp_id == query_snp_id:
                                        continue
                                
                                #Отбор фазированных генотипов по сэмплам.
                                oppos_snp_phased_genotypes = get_phased_genotypes(oppos_snp_data,
                                                                                  sample_names)
                                
                                #Получение значений LD с
                                #помощью оффлайн-калькулятора.
                                #Полезный побочный продукт функции -
                                #частоты альтернативного аллеля
                                #запрашиваемого и кандидатного SNPs.
                                trg_vals = ld_calc(query_snp_phased_genotypes,
                                                   oppos_snp_phased_genotypes)
                                
                                #Отбор кандидатных
                                #SNP по устновленному
                                #исследователем
                                #порогу LD.
                                if trg_vals[thres_ld_measure] < thres_ld_value:
                                        continue
                                
                                #Добавление в конечный список очередного элемента.
                                #Что это будет за элемент - зависит от решения
                                #исследователя, в каком виде выводить результаты.
                                #Для начала проверим, не предпочёл ли исследователь
                                #самый минималистичный вывод - набор refSNPID.
                                #Если да, заселим в список текущий кандидатный
                                #refSNPID и сразу перейдём к словарю со следующим.
                                if trg_file_type == 'rsids':
                                        linked_snps.append(oppos_snp_id)
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
                                oppos_snp_pos = oppos_snp_data['POS']
                                oppos_snp_ref = oppos_snp_data['REF']
                                oppos_snp_alt = oppos_snp_data['ALT']
                                oppos_snp_type = oppos_snp_data['INFO']['VT']
                                oppos_snp_vals = [oppos_snp_pos,
                                                  oppos_snp_id,
                                                  oppos_snp_ref,
                                                  oppos_snp_alt,
                                                  oppos_snp_type,
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
                                        linked_snps.append(dict(zip(oppos_snp_keys,
                                                                    oppos_snp_vals)))
                                        
                                #Если output - TSV - упомянутые значения
                                #сконкатенируются в табличную строку,
                                #которая поступит в конечный список.
                                elif trg_file_type == 'tsv':
                                        linked_snps.append('\t'.join(map(str, oppos_snp_vals)))
                                        
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
                                                                    'gends'], [chrom,
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
                                        query_snp_ref = query_snp_data['REF']
                                        query_snp_alt = query_snp_data['ALT']
                                        query_snp_type = query_snp_data['INFO']['VT']
                                        query_ann_vals, header_keys = [query_snp_pos,
                                                                       query_snp_id,
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
                                                continue
                                        
                                        #TSV: Собираем, во-первых, табулированную шапку
                                        #из названий сниповых характеристик, добавив
                                        #к ним элементы r2, D' и dist, а, во-вторых,
                                        #строку с характеристиками запрашиваемого SNP.
                                        elif trg_file_type == 'tsv':
                                                linked_snps.insert(1, '#' + '\t'.join(header_keys +
                                                                                      ["r2", "D'", 'dist']))
                                                linked_snps.insert(2, '\t'.join(list(map(str, query_ann_vals)) +
                                                                                ['quer'] * 3))
                                                
                                        #RefSNPIDs: дополняем конечный
                                        #список именем снипового столбца и
                                        #идентификатором запрашиваемого SNP.
                                        elif trg_file_type == 'rsids':
                                                linked_snps.insert(1, '#rsID')
                                                linked_snps.insert(2, query_snp_id)
                                                
                                        #Независимо от того, выбран табличный
                                        #формат или вывод одиночных refSNPID,
                                        #прописываем все собранные в конечный
                                        #список данные как обычные строки.
                                        for line in linked_snps:
                                                trg_file_opened.write(line + '\n')
                                                
                        elif verbose in ['yes', 'y']:
                                print(f'\t{thres_ld_measure} >= {thres_ld_value}: SNPs в LD не найдены')
                                
client.close()
