__version__ = 'V4.3'

print('''
Программа ищет в пределах фланков SNPs,
обладающие надпороговым значением неравновесия
по сцеплению с каждым запрашиваемым SNP.

Автор: Платон Быкадоров (platon.work@gmail.com), 2018-2019.
Версия: V4.3.
Лицензия: GNU General Public License version 3.
Поддержать проект: https://money.yandex.ru/to/41001832285976
Документация: https://github.com/PlatonB/ld-tools/blob/master/README.md

Обязательно! Установка модуля:
sudo pip3 install pysam

Поддерживаемые исходные файлы - таблицы,
содержащие столбец с набором refSNPIDs.
Если таких столбцов - несколько,
программа будет использовать самый левый.

Если настройки, запрашиваемые в рамках интерактивного
диалога, вам непонятны - пишите, пожалуйста, в Issues.
''')

def check_input(var):
        '''
        Проверка, правильно ли пользователь ответил на
        запрос, требующий подтверждения или отрицания.
        В случае ошибки работа программы завершится.
        '''
        if var != 'yes' and var != 'y' and var != 'no' \
           and var != 'n' and var != '':
                print(f'{var} - недопустимая опция')
                sys.exit()
                
####################################################################################################

print('\nИмпорт модулей программы...')

import sys

#Подавление формирования питоновского кэша с
#целью предотвращения искажения результатов.
sys.dont_write_bytecode = True

import random, os, re, gzip, dbm, json
sys.path.insert(0, os.path.join(os.getcwd(), 'backend'))
from prepare_intgen_data import process_intgen_data
from retrieve_sample_indices import retrieve_sample_indices
from ld_calc import ld_calc
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
(в любом регистре)
(несколько - через пробел)
(игнорирование ввода ==> всех)
[ALL(|<enter>)|EUR|JPT|JPT AMR YRI|...]: ''').upper().split()
if populations == []:
        populations = ['ALL']
        
genders = input('''\nДля индивидов какого пола считать LD?
(игнорирование ввода ==> обоих)
[male(|m)|female(|f)|both(|<enter>)]: ''').split()
if genders == ['m']:
        genders = ['male']
elif genders == ['f']:
        genders = ['female']
elif genders == [] or genders == ['both']:
        genders = ['male', 'female']
elif genders != ['male'] and genders != ['female']:
        print(f'{genders[0]} - недопустимая опция')
        sys.exit()
        
flank_size = input('''\nРазмер каждого из фланков, в пределах
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
elif thres_ld_measure != 'r_square' and thres_ld_measure != 'd_prime':
        print(f'{thres_ld_measure} - недопустимая опция')
        sys.exit()
        
thres_ld_value = float(input(f'\n{thres_ld_measure} ≥ '))
        
#Вызов функции, которая скачает заархивированные
#VCF и панель сэмплов проекта 1000 Genomes,
#если они ещё не скачаны, а также разместит
#строки изо всех VCF в похромосомные dbm-базы.
#Каждая dbm-база данных (далее - база) будет
#состоять из ключей - refSNPID определённой хромосомы,
#и значений - сжатых строк упомянутых VCF.
#Функция не только скачивает данные, но и
#возвращает абсолютные пути: без расширений -
#к vcf.gz-архивам (они же - к dbm-базам),
#а с расширением - к панели сэмплов.
intgen_sampjson_path, intgen_vcfbase_paths = process_intgen_data(intgen_dir_path)

#Вызов функции, производящей отбор
#сэмплов, относящихся к указанным
#пользователем популяциям и полам,
#и возвращающей индексы этих сэмплов.
sample_indices = retrieve_sample_indices(intgen_sampjson_path,
                                         populations,
                                         genders,
                                         random.choice(intgen_vcfbase_paths) + '.dbm')

#Работа с исходными файлами.
src_file_names = os.listdir(src_dir_path)
for src_file_name in src_file_names:
        src_file_base = '.'.join(src_file_name.split('.')[:-1])
        
        #Построение пути к папке для результатов
        #по текущему пользовательскому файлу,
        #которые, в свою очередь, должны будут
        #распределяться по хромосомным подпапкам.
        trg_dir_path = os.path.join(trg_top_dir_path, f'{src_file_base}_lnkd')
        
        #Открытие файла пользователя на чтение.
        with open(os.path.join(src_dir_path, src_file_name)) as src_file_opened:
                
                print(f'\n\n{src_file_name}')
                
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
                        
                        print(f'\n{query_rs_id}...')
                        
                        #Список, в который будут помещаться,
                        #как минимум, идентификаторы SNP,
                        #сцепленных с запрашиваемым.
                        #Туда также могут попадать
                        #различные характеристики SNPs.
                        #Содержимое списка будет зависить
                        #от пресета, выбранного пользователем.
                        linked_snps = []
                        
                        #Перебор путей к 1000 Genomes-файлам без расширений.
                        for intgen_vcfbase_path in intgen_vcfbase_paths:
                                
                                #К этим путям добавляется
                                #расширение dbm, чтобы можно
                                #было обращаться к базам.
                                #Открытие каждой базы на чтение.
                                with dbm.open(intgen_vcfbase_path + '.dbm') as intgen_vcfdb_opened:
                                        
                                        #Попытка извлечения из очередной базы сжатой строки,
                                        #соответствующей текущему пользовательскому refSNPID.
                                        #Декомпрессия результата, превращение из
                                        #байт-строки в обычную и преобразование в список.
                                        #Если refSNPID в базе нет, то вместо сжатой строки вернётся None.
                                        #В результате разархивации None выйдет пустая байтовая строка.
                                        #После её конвертаций, упомянутых выше, сформируется список
                                        #с обычной пустой строкой в качестве единственного элемента.
                                        query_snp_row = gzip.decompress(intgen_vcfdb_opened.get(query_rs_id)).decode('utf-8').split('\t')
                                        
                                        #Если содержимое списка получилось
                                        #иное, нежели пустая строка, значит
                                        #refSNPID в одной из баз нашёлся.
                                        #После работы с соответствующим
                                        #SNP поиск будет прерван.
                                        if query_snp_row != ['']:
                                                
                                                #Получение номера хромосомы и позиции запрашиваемого SNP.
                                                chr_num, query_snp_pos = query_snp_row[0], int(query_snp_row[1])
                                                
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
                                                trg_file_name = f'{chr_num}_{query_rs_id}_{thres_ld_measure[0]}_{str(thres_ld_value)}.json'
                                                trg_file_path = os.path.join(trg_chrdir_path, trg_file_name)
                                                
                                                #Если запрашиваемый SNP встретился
                                                #повторно, то необходимо предотвратить
                                                #поиск SNPs, с ним сцепленных.
                                                if os.path.exists(trg_file_path) == True:
                                                        print('\tуже был ранее обработан')
                                                        break
                                                
                                                #Открытие на чтение того из tabix-индексированных
                                                #1000 Genomes-VCF, в котором точно есть запрашиваемый SNP.
                                                #Pysam требует на вход лишь путь к vcf.gz-файлу,
                                                #а tbi-сателлита обнаруживает самостоятельно.
                                                #Из этого vcf.gz-файла будут вытаскиваться сцепленные SNPs.
                                                with VariantFile(intgen_vcfbase_path + '.vcf.gz') as variant_file_obj:
                                                        
                                                        #Перебор объектов, соответствующих
                                                        #заданному пользователем координатному
                                                        #окну вокруг запрашиваемого SNP.
                                                        #Каждый объект - это представленная
                                                        #pysam'ом refSNPID-содержащая строка.
                                                        #SNPs из описываемых объектов
                                                        #будем называть кандидатными.
                                                        for rec in variant_file_obj.fetch(chr_num,
                                                                                          query_snp_pos - flank_size,
                                                                                          query_snp_pos + flank_size):
                                                                
                                                                #Преобразование pysam-объекта
                                                                #в строку, а её - в список.
                                                                oppos_snp_row = str(rec).split('\n')[0].split('\t')
                                                                
                                                                #Получение из этого списка таких важных
                                                                #элементов, как позиция, идентификатор
                                                                #и info-ячейка кандидатного SNP.
                                                                #В последней из перечисленных может
                                                                #встречаться флаг MULTI_ALLELIC, позволяющий
                                                                #идентифицировать не-биаллельные SNP.
                                                                oppos_snp_pos, oppos_rs_id, oppos_snp_info = int(oppos_snp_row[1]), \
                                                                                                             oppos_snp_row[2], \
                                                                                                             oppos_snp_row[7]
                                                                
                                                                #Кандидатный SNP имеет шанс войти в конечный
                                                                #список при соответствии нескольким критериям:
                                                                #его идентификатор - refSNPID, он - не мультиаллельный
                                                                #и не тот же самый, что и запрашиваемый.
                                                                if re.match(r'rs\d+$', oppos_rs_id) != None \
                                                                   and oppos_snp_info.find('MULTI_ALLELIC') == -1 \
                                                                   and oppos_rs_id != query_rs_id:
                                                                        
                                                                        #Вычисление LD.
                                                                        ld_vals = ld_calc(query_snp_row,
                                                                                          oppos_snp_row,
                                                                                          sample_indices)
                                                                        
                                                                        #Отбор кандидатных SNP по
                                                                        #пользовательскому порогу LD.
                                                                        if ld_vals[thres_ld_measure] < thres_ld_value:
                                                                                continue
                                                                        
                                                                        #Добавление в конечный список словаря
                                                                        #с позицией и ID сцепленного SNP, а также
                                                                        #со значениями LD и физическим расстоянием
                                                                        #между запрашиваемым и сцепленным SNP.
                                                                        linked_snps.append({'linked.hg38': oppos_snp_pos,
                                                                                            'linked.rsID': oppos_rs_id,
                                                                                            "r2": ld_vals['r_square'],
                                                                                            "D'": ld_vals['d_prime'],
                                                                                            'distance': oppos_snp_pos - query_snp_pos})
                                                                        
                                                #Конечная папка, хромосомная подпапка и
                                                #файл с результатами могут быть созданы
                                                #при условии, что в конечном списке
                                                #оказался хоть один сцепленный SNP.
                                                if linked_snps != []:

                                                        #Если файлу с результатами
                                                        #быть, конечный список дополнится
                                                        #первым элементом, содержащим
                                                        #номер хромосомы, позицию
                                                        #и ID запрашиваемого SNP,
                                                        #а также параметры запроса.
                                                        linked_snps.insert(0, {'chr': int(chr_num),
                                                                               'queried.hg38': query_snp_pos,
                                                                               'queried.rsID': query_rs_id,
                                                                               'flank': flank_size,
                                                                               thres_ld_measure: thres_ld_value,
                                                                               'populations': populations,
                                                                               'genders': genders})
                                                        
                                                        #Создание конечной папки и хромосомной
                                                        #подпапки, если таковых ещё нет.
                                                        #Если подпапка есть, то наличие
                                                        #папки проверяться не будет.
                                                        if os.path.exists(trg_chrdir_path) == False:
                                                                if os.path.exists(trg_dir_path) == False:
                                                                        os.mkdir(trg_dir_path)
                                                                os.mkdir(trg_chrdir_path)

                                                        #Создание и открытие конечного файла на запись.
                                                        with open(trg_file_path, 'w') as trg_file_opened:
                                                                
                                                                #Полученный ранее список словарей пропишется
                                                                #в JSON-файл с формированием отступов.
                                                                json.dump(linked_snps, trg_file_opened, indent=4)
                                                                
                                                else:
                                                        print(f'\t{thres_ld_measure} >= {thres_ld_value}: SNPs в LD не найдены')
                                                        
                                                break
                                        
                        #Если refSNPID из файла
                        #пользователя не обнаружился
                        #в базах, значит он - невалидный.
                        #В конечный файл он допущен не будет.
                        #Одна из вероятных причин - в том,
                        #что именуемый им SNP - не биаллельный.
                        else:
                                print('\tневалидный refSNPID (возможно, ID мультиаллельного SNP).')
                                continue
